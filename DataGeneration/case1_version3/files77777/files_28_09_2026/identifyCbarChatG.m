function R = identifyCbarChatG(S, cfg)
% IDENTIFYCBARCHATG  Identify Cbar, Chat and G from energy, integrated stress
%                    and integrated double stress, one set for every N.
%
%   R = identifyCbarChatG(S, cfg)     S from smoothData.m
%
% ENERGY DENSITY (C = C_hom known, psi = sym(psi), gamma = eps_M - psi,
% kappa_ijk = d psi_ij / d x_k, ell = L/N; Voigt [a11 a22 2a12]):
%
%   W = 1/2 psi'C psi + 1/2 gamma'Cbar gamma + psi'Chat gamma
%       + 1/2 ell^2 kappa' diag(G11,G22,G22,G11,G33,G33) kappa
%
%   Cbar = [Cb11 Cb12 0; Cb12 Cb11 0; 0 0 Cb33]   (Cbar_11 = Cbar_22)
%   Chat = [Ch11 Ch12 0; Ch12 Ch11 0; 0 0 Ch33]   (Chat_11 = Chat_22)
%   cfg.id.Gmode = 'single': G11 = G22 = G33 = G  -> 7 unknowns
%                  'three' : G11, G22, G33         -> 9 unknowns
%   kappa order [k111 k112 k221 k222 k121 k122]
%
% Conjugate stresses (Voigt stress [s11 s22 s12]):
%   s_(sym)       = C psi + Chat gamma
%   sigma_(gamma) = Cbar gamma + Chat psi        (conjugate to eps_M)
%   mu            = ell^2 diag(G..) kappa
%
% EQUATIONS for every load case at the given N (A_i, I_i as in dataGeneration):
%  energy (1)   C11 A1 + 2C12 A2 + C22 A3 + 4C33 A4
%               + Cb11 (A5+A7) + 2 Cb12 A6 + 4 Cb33 A8
%               + 2 Ch11 (A9+A12) + 2 Ch12 (A10+A11) + 8 Ch33 A13
%               + ell^2 G (A14+...+A19)                     = 2U_het
%  stress (3)   int sigma_(gamma):  [Cb11 I4 + Cb12 I5 + Ch11 I1 + Ch12 I2;
%                                    Cb12 I4 + Cb11 I5 + Ch12 I1 + Ch11 I2;
%                                    2 Cb33 I6 + 2 Ch33 I3]  = int sigma_FE
%  double (6)   c_ij int sigma_(gamma),ij x_k + int mu_ijk  = c_ij int sigma_FE,ij x_k
%               with c = 1 for 11, 22 and c = 2 for 12 (sigma_12 and sigma_21),
%               int mu = ell^2 G [I7..I12]
%               cfg.id.doubleStressModel = 'mu' keeps only int mu on the left.
%
% Each group is divided by its size (2U_het, |int sigma|, |moment vector|),
% weighted by wE / wS / wD, and used only if informative (it vanishes by
% symmetry for pure G or pure H loads).
%
% MODE 'anchored' (default):
%   Cbar(N) = C_inf + sum_k Dbar_k N^-e_k,  Chat(N) = C_inf + sum_k Dhat_k N^-e_k,
%   C_inf = C_hom (both tend to the homogenized stiffness for N -> inf),
%   e = cfg.id.CExponents (default [1 2]), deviation signs free or restricted
%   (cfg.id.CbDeviation / ChDeviation = 'free' | 'decreasing' | 'increasing');
%   G11, G22, G33 CONSTANT (independent of N, cfg.id.Gmode = 'three'), >= 0;
%   curvature energy 1/2 ell^m kappa' diag(G11,G22,G22,G11,G33,G33) kappa,
%   ell = L/N, m fixed (cfg.id.m = 2) or identified (cfg.id.fitM = true: the
%   convex problem is solved for every trial m, grid + golden section).
%   Constraints at every N: Cbar > 0, Chat > 0 (margin), energy block
%   chat_i^2 <= (1 - blockMargin) c_i cbar_i (semi-definite: in the limit
%   Cbar = Chat = C_hom the local energy is the classical 1/2 eps_M'C eps_M).
%
% MODE 'monotone': every parameter is a monotone function of N,
%   p_j(N) = sum_k c_jk N^-e_k,  c_jk >= 0,  e = cfg.id.trendExponents (0 1 2)
% so all components are positive and non-increasing in N. All N are fitted
% together, subject at every N to
%   Chat positive definite:        chat_i >= m                       (linear)
%   whole local energy positive definite, [C Chat; Chat Cbar] > 0:
%       cbar_i * c_i >= chat_i^2 / (1 - delta)                       (convex)
% where c_i, chat_i, cbar_i are the eigenvalues of C, Chat, Cbar in the common
% eigenvectors (1,1,0)/sqrt2, (1,-1,0)/sqrt2, (0,0,1) of the cubic structure:
%   mode 1: p11 + p12,  mode 2: p11 - p12,  mode 3: p33.
% The problem is convex (quadratic objective, convex quadratic constraints);
% it is solved with fmincon (sqp) or, in Octave, with sqp.
%
% MODES 'perN' / 'global' (free in N; kept for comparison):
% CONSTRAINTS (linear, because of the cubic structure):
%   Cbar pos. def.:  Cb11 - Cb12 >= m, Cb11 + Cb12 >= m, Cb33 >= m
%   Chat pos. def.:  Ch11 - Ch12 >= m, Ch11 + Ch12 >= m, Ch33 >= m
%   G >= 0,     m = pdMargin * max eig(C)
% The convex quadratic program is solved exactly by enumerating active sets.
% Optionally (requireEnergyPD) the whole local energy block [C Chat; Chat Cbar]
% is kept positive definite as well (non-linear constraint, fmincon / sqp).

id = cfg.id;
C  = S.C_hom;
C  = [C(1,1) C(1,2) 0; C(1,2) C(2,2) 0; 0 0 C(3,3)];     % C11 C12 C22 C33 as in W
lamC = max(eig(C));
nC = size(S.A, 2);
cases = id.cases;  if isempty(cases), cases = 1:nC; end
Nk = S.Nlist(:);  nK = numel(Nk);
switch lower(id.Gmode)
    case 'single', pNames = {'Cb11','Cb12','Cb33','Ch11','Ch12','Ch33','G'};
    case 'three',  pNames = {'Cb11','Cb12','Cb33','Ch11','Ch12','Ch33','G11','G22','G33'};
    otherwise, error('cfg.id.Gmode must be ''single'' or ''three''.');
end
nPar = numel(pNames);
[Gc, h, cNames] = pdConstraints(nPar, id, lamC, pNames);
chName = {'energy', 'stress', 'double stress'};

R = struct('paramNames', {pNames}, 'channelNames', {chName}, 'Nlist', Nk.', 'cases', cases, ...
           'C', C, 'Gmode', id.Gmode, 'doubleStressModel', id.doubleStressModel, ...
           'weights', [id.wE id.wS id.wD]);

if ~isfield(id, 'trendExponents'), id.trendExponents = [0 1 2]; end
if ~isfield(id, 'energyPD'),       id.energyPD = true;         end
if ~isfield(id, 'compareFree'),    id.compareFree = true;      end
if ~isfield(id, 'm'),              id.m = 2;                   end
if ~isfield(id, 'fitM'),           id.fitM = false;            end
if ~isfield(id, 'mRange'),         id.mRange = [0.5 4];        end
if ~isfield(id, 'mGridN'),         id.mGridN = 15;             end
if ~isfield(id, 'CExponents'),     id.CExponents = [1 2];      end
if ~isfield(id, 'CbDeviation'),    id.CbDeviation = 'free';    end
if ~isfield(id, 'ChDeviation'),    id.ChDeviation = 'free';    end
if ~isfield(id, 'blockMargin'),    id.blockMargin = 0;         end
doMono = strcmpi(id.mode, 'monotone');
doAnch = strcmpi(id.mode, 'anchored');

% ---------------- anchored: C's -> C_hom, G constant, ell^m ----------------
if doAnch
    M = solveAnchoredProfile(S, cases, id, cfg, C, pNames);
    R.anchored = M;
    id.m = M.m;                       % comparison runs use the same length-scale exponent
    R.m  = M.m;
    printAnchored(M, C, lamC, pNames);
end

% ---------------- monotone in N, all N together ----------------
if doMono
    M = solveMonotone(S, cases, id, cfg, C, pNames);
    R.monotone = M;
    printMonotone(M, C, lamC, pNames);
end

% ---------------- per N (free in N) ----------------
if any(strcmpi(id.mode, {'perN', 'both'})) || ((doMono || doAnch) && id.compareFree)
    if doMono || doAnch
        fprintf('\n--- for comparison: each N identified on its own (C''s and G free per N, m = %g, energy block PD not imposed) ---\n', id.m);
    end
    P = nan(nK, nPar);  rms = nan(nK, 3);  info = zeros(nK, 3);  cnd = nan(nK,1);
    act = cell(nK,1);  meth = cell(nK,1);  eigs3 = nan(nK, 3);
    for iN = 1:nK
        r = solveBlock(S, cases, iN, id, cfg, C, Gc, h, cNames);
        P(iN,:) = r.theta;  rms(iN,:) = r.rms;  info(iN,:) = r.nInfo;  cnd(iN) = r.cond;
        act{iN} = r.active;  meth{iN} = r.method;  eigs3(iN,:) = r.minEig;
        if ~isempty(r.unseen)
            fprintf('  N = %d: no information on %s (set by the bounds only)\n', Nk(iN), strjoin(pNames(r.unseen), ', '));
        end
    end
    R.perN = struct('params', P, 'rmsPerChannel', rms, 'nInformative', info, 'cond', cnd, ...
                    'activeConstraints', {act}, 'method', {meth}, 'minEig', eigs3);
    R.trend = fitTrends(Nk, P, pNames);
    printPerN(R, C, lamC);
end

% ---------------- one set for all N ----------------
if any(strcmpi(id.mode, {'global', 'both'}))
    r = solveBlock(S, cases, 1:nK, id, cfg, C, Gc, h, cNames);
    R.global = struct('params', r.theta, 'rmsPerChannel', r.rms, 'activeConstraints', {r.active}, ...
                      'method', r.method, 'minEig', r.minEig, 'cond', r.cond);
    fprintf('\n[identifyCbarChatG] one set for all N = %s (%s):\n', mat2str(Nk.'), r.method);
    for j = 1:nPar, fprintf('    %-5s = %12.6g\n', pNames{j}, r.theta(j)); end
    fprintf('  RMS misfit: energy %.4f | stress %.4f | double stress %.4f\n', r.rms);
end

% ---------------- output ----------------
if doAnch
    R.params = R.anchored.params;
elseif doMono
    R.params = R.monotone.params;  R.trendCoef = R.monotone.coef;
elseif isfield(R, 'perN')
    R.params = R.perN.params;
end
save(fullfile(cfg.outDir, 'IdentifiedCbarChatG.mat'), '-struct', 'R');
writeCSV(fullfile(cfg.outDir, 'IdentifiedCbarChatG_vsN.csv'), R);
if cfg.plots && doAnch
    plotAnchored(R, C, lamC);
    plotFit(S, R.anchored.params, cases, id, cfg, C, 'Fit_anchored', ...
        sprintf('Cbar(N), Chat(N) -> C_{hom}, constant G, m = %.3g', R.anchored.m));
elseif cfg.plots && doMono
    plotMonotone(R, C, lamC);
    plotFit(S, R.monotone.params, cases, id, cfg, C, 'Fit_monotone', ...
        'monotone Cbar(N), Chat(N), G(N), energy positive definite');
elseif cfg.plots
    if isfield(R, 'perN')
        plotParams(R, C, lamC);
        plotFit(S, R.perN.params, cases, id, cfg, C, 'Fit_perN', 'Cbar, Chat, G identified at each N');
    end
    if isfield(R, 'global')
        plotFit(S, repmat(R.global.params, nK, 1), cases, id, cfg, C, 'Fit_global', ...
            'one Cbar, Chat, G for all N');
    end
end
end


% ============================================================
% ============== equations and constrained solve =============
% ============================================================

function [Prow, brow, chan] = caseRows(S, ic, iN, id, cfg, C)
    % weighted, scaled rows of one load case at one N
    V  = S.VolN(iN);                 % volume of the region (full or interior) at this N
    Lr = S.LrN(iN);                  % rms moment arm of the region
    single = strcmpi(id.Gmode, 'single');
    mExp = 2;  if isfield(id, 'm'), mExp = id.m; end
    a  = squeeze(S.A(iN,ic,:));  I = squeeze(S.I(iN,ic,:));
    l2 = S.ell2(iN)^(mExp/2);        % ell^m  (ell^2 for m = 2)
    U2 = S.twoU_het(iN,ic);
    sigRef = sqrt(U2*max(eig(C))/V) * V;
    Prow = [];  brow = [];  chan = [];

    % ---- energy ----
    known = C(1,1)*a(1) + 2*C(1,2)*a(2) + C(2,2)*a(3) + 4*C(3,3)*a(4);
    ce = [a(5)+a(7), 2*a(6), 4*a(8), 2*(a(9)+a(12)), 2*(a(10)+a(11)), 8*a(13)];
    if single
        cg = l2*sum(a(14:19));
    else   % A = diag(G11,G22,G22,G11,G33,G33): G11 ~ k111,k222; G22 ~ k112,k221; G33 ~ k121,k122
        cg = l2*[a(14)+a(17), a(15)+a(19), a(18)+a(16)];
    end
    w = sqrt(id.wE) / U2;
    Prow = [Prow; w*[ce, cg]];  brow = [brow; w*(U2 - known)];  chan = [chan; 1];

    % ---- integrated stress ----
    Sfe = squeeze(S.Sfe(iN,ic,:));
    nG  = 1 + 2*(~single);
    if norm(Sfe) >= id.tolInfo*sigRef
        blk = [ I(4) I(5) 0      I(1) I(2) 0;
                I(5) I(4) 0      I(2) I(1) 0;
                0    0    2*I(6) 0    0    2*I(3) ];
        w = sqrt(id.wS) / norm(Sfe);
        Prow = [Prow; w*[blk, zeros(3,nG)]];  brow = [brow; w*Sfe];  chan = [chan; 2*ones(3,1)];
    end

    % ---- integrated double stress ----
    Mfe = squeeze(S.Mfe(iN,ic,:));
    m   = Mfe .* [1;1;1;1;2;2];                         % conjugate to kappa
    if norm(m) >= id.tolInfo*sigRef*Lr
        loc = [ I(19) I(21) 0       I(13) I(15) 0;
                I(20) I(22) 0       I(14) I(16) 0;
                I(21) I(19) 0       I(15) I(13) 0;
                I(22) I(20) 0       I(16) I(14) 0;
                0     0     4*I(23) 0     0     4*I(17);
                0     0     4*I(24) 0     0     4*I(18) ];
        if strcmpi(id.doubleStressModel, 'mu'), loc(:) = 0; end
        k = I(7:12);
        if single
            gm = l2*k;
        else
            gm = l2*[ [k(1);0;0;k(4);0;0], [0;k(2);k(3);0;0;0], [0;0;0;0;k(5);k(6)] ];
        end
        w = sqrt(id.wD) / norm(m);
        Prow = [Prow; w*[loc, gm]];  brow = [brow; w*m];  chan = [chan; 3*ones(6,1)];
    end
end

function r = solveBlock(S, cases, iNs, id, cfg, C, Gc, h, cNames)
    P = [];  b = [];  chan = [];  nInfo = zeros(1,3);
    for ic = cases(:).'
        for iN = iNs(:).'
            [Pr, br, ch] = caseRows(S, ic, iN, id, cfg, C);
            P = [P; Pr];  b = [b; br];  chan = [chan; ch];               %#ok<AGROW>
            nInfo = nInfo + [1, any(ch == 2), any(ch == 3)];
        end
    end
    wq = [id.wE id.wS id.wD];
    use = wq(chan).' > 0;
    Pu = P(use,:);  bu = b(use);
    s = max(abs(Pu), [], 1);                                          % column scaling, floored so
    s = 1 ./ max(s, 1e-6*max(s));                                     % unseen parameters stay well posed
    Ps = Pu .* s;  Gs = Gc .* s;
    [x, active, method] = qpActiveSet(Ps, bu, Gs, h, cNames);
    theta = (x(:) .* s(:)).';
    theta(7:end) = max(theta(7:end), 0);                              % G >= 0 (round-off)
    if id.requireEnergyPD
        [theta, method, active] = energyPDsolve(theta, Pu, bu, Gc, h, C, id);
    end
    res = P*theta(:) - b;
    rms = nan(1,3);
    for q = 1:3
        mq = use & chan == q;
        if any(mq), rms(q) = sqrt(mean((res(mq)/sqrt(wq(q))).^2)); end
    end
    sv = svd(Ps);
    colMax = max(abs(Pu), [], 1);
    unseen = find(colMax < 1e-6*max(colMax));                         % data carry no information
    r = struct('theta', theta, 'rms', rms, 'nInfo', nInfo, 'cond', sv(1)/max(sv(end), eps), ...
               'active', {active}, 'method', method, 'minEig', minEigs(theta, C), 'unseen', unseen);
end

function M = solveAnchoredProfile(S, cases, id, cfg, C, pNames)
    % m fixed (default 2), or identified by profiling the convex problem over m
    if ~id.fitM
        M = solveAnchored(S, cases, id, cfg, C, pNames, id.m);
        M.mProfile = [id.m, M.objective];  M.mFitted = false;  M.mFlat = false;
        return;
    end
    Jm = @(mm) getfield(solveAnchored(S, cases, id, cfg, C, pNames, mm), 'objective');
    mg = linspace(id.mRange(1), id.mRange(2), id.mGridN);
    Jg = arrayfun(Jm, mg);
    prof = [mg(:), Jg(:)];
    [~, q0] = min(Jg);
    a = mg(max(q0-1, 1));  b = mg(min(q0+1, numel(mg)));
    gr = (sqrt(5) - 1)/2;                               % golden-section refinement on [a, b]
    c1 = b - gr*(b - a);  c2 = a + gr*(b - a);  f1 = Jm(c1);  f2 = Jm(c2);
    prof = [prof; c1 f1; c2 f2];
    for it = 1:25
        if f1 < f2
            b = c2;  c2 = c1;  f2 = f1;  c1 = b - gr*(b - a);  f1 = Jm(c1);  prof(end+1,:) = [c1 f1]; %#ok<AGROW>
        else
            a = c1;  c1 = c2;  f1 = f2;  c2 = a + gr*(b - a);  f2 = Jm(c2);  prof(end+1,:) = [c2 f2]; %#ok<AGROW>
        end
        if b - a < 1e-3, break; end
    end
    M = solveAnchored(S, cases, id, cfg, C, pNames, 0.5*(a + b));
    M.mProfile = sortrows(prof);  M.mFitted = true;
    M.mFlat = (max(Jg) - min(Jg)) / max(min(Jg), realmin) < 1e-3;
end

function M = solveAnchored(S, cases, id, cfg, C, pNames, m)
    % All N together:  X(N) = X_inf + sum_k D_k N^-e_k  for X = Cbar, Chat (X_inf = C_hom),
    % G constant, curvature energy 1/2 ell^m kappa' diag(G..) kappa, ell = L/N.
    id.m = m;
    Nk = S.Nlist(:);  nK = numel(Nk);  np = numel(pNames);  nGp = np - 6;
    e  = id.CExponents(:).';  nb = numel(e);  nx = 6*nb + nGp;
    aC = 0.5*(C(1,1) + C(2,2));
    th0 = [aC, C(1,2), C(3,3), aC, C(1,2), C(3,3), zeros(1, nGp)];      % limit N -> inf
    Pb = [];  bb = [];  chan = [];  nIdx = [];  Tn = cell(nK, 1);
    for iN = 1:nK
        T = zeros(np, nx);
        T(1:6, 1:6*nb) = kron(eye(6), Nk(iN).^(-e));    % deviations of Cbar, Chat
        T(7:np, 6*nb+1:nx) = eye(nGp);                  % G: constant in N
        Tn{iN} = T;
        for ic = cases(:).'
            [Pr, br, ch] = caseRows(S, ic, iN, id, cfg, C);
            Pb = [Pb; Pr*T];  bb = [bb; br - Pr*th0.'];  chan = [chan; ch];            %#ok<AGROW>
            nIdx = [nIdx; iN*ones(numel(br), 1)];                                        %#ok<AGROW>
        end
    end
    wq  = [id.wE id.wS id.wD];  use = wq(chan).' > 0;
    Pu  = Pb(use,:);  bu = bb(use);
    sc  = max(abs(Pu), [], 1);  sc = 1 ./ max(sc, 1e-6*max(sc));          % x = x_s .* sc
    Ps  = Pu .* sc;

    % eigen-modes (common eigenvectors of the cubic matrices)
    cmo = [aC + C(1,2); aC - C(1,2); C(3,3)];
    modeP = [1 1 0; 1 -1 0; 0 0 1];
    mrg = id.pdMargin * max(eig(C));  del = id.blockMargin;
    Uh = zeros(3*nK, nx);  Vb = Uh;  cI = zeros(3*nK, 1);
    for iN = 1:nK
        for i = 1:3
            u = zeros(1, np);  u(4:6) = modeP(i,:);  v = zeros(1, np);  v(1:3) = modeP(i,:);
            r = 3*(iN-1) + i;
            Uh(r,:) = u*Tn{iN};  Vb(r,:) = v*Tn{iN};  cI(r) = cmo(i);
        end
    end
    % mode values: chat_i(N) = c_i + Uh x,  cbar_i(N) = c_i + Vb x
    Alin = [-Uh; -Vb] .* sc;  blin = [cI - mrg; cI - mrg];     % Chat > 0, Cbar > 0 (margin)

    % bounds: sign of the deviations (optional), G >= 0
    lb = -inf(nx,1);  ub = inf(nx,1);
    devIdx = @(jj) reshape((jj(:)-1)*nb + (1:nb), 1, []);
    switch lower(id.CbDeviation)
        case 'decreasing', lb(devIdx(1:3)) = 0;
        case 'increasing', ub(devIdx(1:3)) = 0;
    end
    switch lower(id.ChDeviation)
        case 'decreasing', lb(devIdx(4:6)) = 0;
        case 'increasing', ub(devIdx(4:6)) = 0;
    end
    lb(6*nb+1:nx) = 0;

    % strictly feasible start: Cbar above C_hom (or Chat below), Chat = C_hom, G small
    x0 = zeros(nx, 1);
    if ~strcmpi(id.CbDeviation, 'increasing')
        x0(devIdx(1:3)) = kron(0.5*th0(1:3), [1, zeros(1, nb-1)]);
    elseif ~strcmpi(id.ChDeviation, 'decreasing')
        x0(devIdx(4:6)) = kron(-0.3*th0(4:6), [1, zeros(1, nb-1)]);
    end
    x0(6*nb+1:nx) = 1e-2*max(eig(C));
    xs0 = x0 ./ sc(:);

    obj = @(xs) monoObjective(xs, Ps, bu);
    if id.energyPD
        nonl = @(xs) anchNonlcon(xs, Uh .* sc, Vb .* sc, cI, del);
    else
        nonl = [];
    end
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
            'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
            'MaxIterations', 5000, 'MaxFunctionEvaluations', 1e6, ...
            'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-15);
        [xs, ~, flag] = fmincon(obj, xs0, Alin, blin, [], [], lb./sc(:), ub./sc(:), nonl, opts);
        method = sprintf('fmincon (sqp), exitflag %d', flag);
    elseif exist('sqp', 'file') == 2 || exist('sqp', 'builtin') == 5
        if isempty(nonl)
            hf = {@(z) blin - Alin*z, @(z) -Alin};
        else
            hf = {@(z) [blin - Alin*z; -firstOut(nonl, z)], @(z) [-Alin; -thirdOutT(nonl, z)]};
        end
        [xs, ~, info] = sqp(xs0, {@(z) firstOut(obj, z), @(z) secondOut(obj, z)}, [], hf, ...
                            lb./sc(:), ub./sc(:), 3000, 1e-12);
        method = sprintf('sqp (Octave), info %d', info);
    else
        error('solveAnchored needs fmincon (Optimization Toolbox) or Octave''s sqp.');
    end

    x   = sc(:) .* xs(:);
    Dev = reshape(x(1:6*nb), nb, 6);                   % Dev(k, j): parameter j, exponent e(k)
    Gv  = x(6*nb+1:nx).';
    P   = zeros(nK, np);
    for iN = 1:nK, P(iN,:) = th0 + (Tn{iN}*x).'; end
    res = Pb*x - bb;
    rms = nan(nK, 3);
    for iN = 1:nK
        for q = 1:3
            mq = use & chan == q & nIdx == iN;
            if any(mq), rms(iN,q) = sqrt(mean((res(mq)/sqrt(wq(q))).^2)); end
        end
    end
    mineig = zeros(nK, 3);
    for iN = 1:nK, mineig(iN,:) = minEigs(P(iN,:), C); end
    M = struct('params', P, 'Dev', Dev, 'G', Gv, 'm', m, 'exponents', e, 'Nlist', Nk.', ...
               'anchor', th0, 'rmsPerChannel', rms, 'minEig', mineig, 'method', method, ...
               'objective', 0.5*norm(res(use))^2);
    M.eval = @(Nq) [repmat(th0(1:6), numel(Nq), 1) + (Nq(:).^(-e))*Dev, repmat(Gv, numel(Nq), 1)];
end

function [c, ceq, gc, gceq] = anchNonlcon(x, U, V, cI, del)
    % energy block positive semi-definite in every mode and N:
    %   chat_i^2 - (1-del) c_i cbar_i <= 0   (convex; scaled by c_i^2)
    ch = cI + U*x(:);  cb = cI + V*x(:);
    c   = (ch.^2 - (1 - del)*cI.*cb) ./ cI.^2;
    ceq = [];
    if nargout > 2
        gc   = ((2*ch./cI.^2) .* U - ((1 - del)./cI) .* V).';
        gceq = [];
    end
end

function printAnchored(M, C, lamC, pN)
    Nk = M.Nlist(:);  P = M.params;
    fprintf('\n[identifyCbarChatG] anchored: Cbar(N), Chat(N) -> C_hom for N -> inf, G constant (%s)\n', M.method);
    if M.mFitted
        fprintf('  length-scale exponent m identified: m = %.3f  (energy 1/2 ell^m kappa''G kappa)\n', M.m);
        if M.mFlat, fprintf('  WARNING: objective almost flat in m -- m is not identifiable from these data\n'); end
    else
        fprintf('  length-scale exponent m fixed: m = %g\n', M.m);
    end
    fprintf('  constant curvature moduli:');
    for j = 7:numel(pN), fprintf('  %s = %.6g', pN{j}, M.G(j-6)); end
    fprintf('\n  C_hom: C11 = %.6g, C12 = %.6g, C22 = %.6g, C33 = %.6g\n', C(1,1), C(1,2), C(2,2), C(3,3));
    fprintf('  %4s', 'N');  fprintf(' %11s', pN{1:6});
    fprintf(' | %8s %8s %8s | %7s %7s %7s\n', 'eigCb', 'eigCh', 'eigTot', 'rmsE', 'rmsS', 'rmsD');
    for iN = 1:numel(Nk)
        fprintf('  %4d', Nk(iN));  fprintf(' %11.5g', P(iN,1:6));
        fprintf(' | %8.3g %8.3g %8.3g | %7.4f %7.4f %7.4f\n', M.minEig(iN,:)/lamC, M.rmsPerChannel(iN,:));
    end
    fprintf('  normalized by C_hom:   N   Cb11/C11  Cb12/C12  Cb33/C33 | Ch11/C11  Ch12/C12  Ch33/C33\n');
    for iN = 1:numel(Nk)
        fprintf('                      %4d %9.4f %9.4f %9.4f | %9.4f %9.4f %9.4f\n', Nk(iN), ...
            P(iN,1)/C(1,1), P(iN,2)/C(1,2), P(iN,3)/C(3,3), P(iN,4)/C(1,1), P(iN,5)/C(1,2), P(iN,6)/C(3,3));
    end
    fprintf('  X(N) = X_inf + sum_k D_k N^-e_k,  e = %s,  X_inf = C_hom:\n', mat2str(M.exponents));
    for j = 1:6
        parts = arrayfun(@(k) termText(M.Dev(k,j), M.exponents(k)), 1:numel(M.exponents), 'UniformOutput', false);
        fprintf('    %-5s(N) = %.6g %s\n', pN{j}, M.anchor(j), strtrim(strjoin(parts, ' ')));
    end
    fprintf('  (eig: smallest eigenvalue of Cbar, Chat, [C Chat; Chat Cbar] / max eig(C); the block is only\n');
    fprintf('   semi-definite in the limit N -> inf, where the energy reduces to the classical one)\n');
end

function plotAnchored(R, C, lamC)
    M = R.anchored;  N = M.Nlist(:);  P = M.params;
    Nf = linspace(min(N), max(N), 300).';  Pf = M.eval(Nf);
    nrm = [C(1,1) C(1,2) C(3,3) C(1,1) C(1,2) C(3,3)];
    free = isfield(R, 'perN');
    figure('Color','w', 'Name', 'Params_vsN_anchored', 'Position', [30 50 1700 900]);
    grp = {1:3, 4:6};  ttl = {'(a) Cbar(N) / C_{hom}', '(b) Chat(N) / C_{hom}'};
    lb  = {{'Cbar_{11}/C_{11}', 'Cbar_{12}/C_{12}', 'Cbar_{33}/C_{33}'}, ...
           {'Chat_{11}/C_{11}', 'Chat_{12}/C_{12}', 'Chat_{33}/C_{33}'}};
    mk = {'o','s','^'};  cl = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
    for g = 1:2
        ax = subplot(2,3,g);  hold(ax,'on');  grid(ax,'on');  hh = [];
        for q = 1:3
            j = grp{g}(q);
            hh(end+1) = plot(ax, N, P(:,j)/nrm(j), mk{q}, 'Color', cl(q,:), 'MarkerFaceColor', cl(q,:)); %#ok<AGROW>
            plot(ax, Nf, Pf(:,j)/nrm(j), '-', 'Color', cl(q,:), 'LineWidth', 1.3);
            if free, plot(ax, N, R.perN.params(:,j)/nrm(j), mk{q}, 'Color', cl(q,:)); end
        end
        plot(ax, [min(N) max(N)], [1 1], 'k:');
        legend(ax, hh, lb{g}, 'Location', 'best');  set(ax, 'XTick', N);  xlabel(ax, 'N');  title(ax, ttl{g});
    end
    ax = subplot(2,3,3);  hold(ax,'on');  grid(ax,'on');  hh = [];  lg = {};
    for j = 7:numel(R.paramNames)
        q = j - 6;
        hh(end+1) = plot(ax, [min(N) max(N)], M.G(q)*[1 1], '-', 'Color', cl(q,:), 'LineWidth', 1.5); %#ok<AGROW>
        if free, plot(ax, N, R.perN.params(:,j), mk{q}, 'Color', cl(q,:)); end
        lg{end+1} = sprintf('%s = %.4g', R.paramNames{j}, M.G(q)); %#ok<AGROW>
    end
    legend(ax, hh, lg, 'Location', 'best');  set(ax, 'XTick', N);  xlabel(ax, 'N');
    title(ax, sprintf('(c) G constant (line); open = each N on its own; m = %.3g', M.m));
    ax = subplot(2,3,4);  hold(ax,'on');  grid(ax,'on');
    plot(ax, N, M.minEig/lamC, '-o');  plot(ax, [min(N) max(N)], [0 0], 'k:');
    legend(ax, {'min eig Cbar', 'min eig Chat', 'min eig [C Chat; Chat Cbar]'}, 'Location', 'best');
    set(ax, 'XTick', N);  xlabel(ax, 'N');  ylabel(ax, '/ max eig(C_{hom})');  title(ax, '(d) positive definiteness');
    ax = subplot(2,3,5);  hold(ax,'on');  grid(ax,'on');
    plot(ax, N, M.rmsPerChannel, '-o');
    legend(ax, {'energy', 'stress', 'double stress'}, 'Location', 'best');
    set(ax, 'XTick', N);  xlabel(ax, 'N');  title(ax, '(e) RMS misfit per channel');
    ax = subplot(2,3,6);  hold(ax,'on');  grid(ax,'on');
    plot(ax, M.mProfile(:,1), M.mProfile(:,2), '.-');  plot(ax, M.m, M.objective, 'ro', 'MarkerFaceColor', 'r');
    xlabel(ax, 'm');  ylabel(ax, 'objective J');  title(ax, '(f) objective versus length-scale exponent m');
    sgtitle('Cbar(N), Chat(N) \rightarrow C_{hom} (N \rightarrow \infty), constant G, curvature energy \propto \ell^m');
end

function M = solveMonotone(S, cases, id, cfg, C, pNames)
    % all N together, p_j(N) = sum_k c_jk N^-e_k, c >= 0, energy PD at every N
    Nk = S.Nlist(:);  nK = numel(Nk);  np = numel(pNames);
    e  = id.trendExponents(:).';  nb = numel(e);  nc = np*nb;
    Bn = Nk .^ (-e);                                   % nK x nb
    Pb = [];  bb = [];  chan = [];  nIdx = [];
    for iN = 1:nK
        TN = kron(eye(np), Bn(iN,:));                  % theta(N) = TN * c  (c param-major)
        for ic = cases(:).'
            [Pr, br, ch] = caseRows(S, ic, iN, id, cfg, C);
            Pb = [Pb; Pr*TN];  bb = [bb; br];  chan = [chan; ch];                 %#ok<AGROW>
            nIdx = [nIdx; iN*ones(numel(br),1)];                                  %#ok<AGROW>
        end
    end
    wq  = [id.wE id.wS id.wD];  use = wq(chan).' > 0;
    Pu  = Pb(use,:);  bu = bb(use);
    sc  = max(abs(Pu), [], 1);  sc = 1 ./ max(sc, 1e-6*max(sc));     % c = sc .* x
    Ps  = Pu .* sc;

    % common eigen-modes of C, Chat, Cbar (cubic structure)
    aC  = 0.5*(C(1,1) + C(2,2));
    cmo = [aC + C(1,2); aC - C(1,2); C(3,3)];          % eigenvalues of C
    modeP = [1 1 0; 1 -1 0; 0 0 1];
    mrg = id.pdMargin * max(eig(C));  del = id.pdMargin;
    Uh = zeros(3*nK, nc);  Vb = zeros(3*nK, nc);  cI = zeros(3*nK, 1);
    for iN = 1:nK
        for i = 1:3
            u = zeros(1,np);  u(4:6) = modeP(i,:);      % Chat mode i
            v = zeros(1,np);  v(1:3) = modeP(i,:);      % Cbar mode i
            r = 3*(iN-1) + i;
            Uh(r,:) = kron(u, Bn(iN,:));  Vb(r,:) = kron(v, Bn(iN,:));  cI(r) = cmo(i);
        end
    end
    Uhx = Uh .* sc;  Vbx = Vb .* sc;                   % in scaled variables x

    % feasible start: Chat = C, Cbar = 2 C, G small (on the lowest exponent)
    th0 = [2*aC, 2*C(1,2), 2*C(3,3), aC, C(1,2), C(3,3), 1e-3*max(eig(C))*ones(1, np-6)];
    th0 = max(th0, 1e-3*max(eig(C)));
    [~, k0] = min(e);  c0 = zeros(nb, np);  c0(k0,:) = th0 * max(Nk)^e(k0);
    x0  = c0(:) ./ sc(:);
    lb  = zeros(nc, 1);  ub = inf(nc, 1);

    obj  = @(x) monoObjective(x, Ps, bu);
    Alin = -Uhx;  blin = -mrg*ones(3*nK,1);            % Chat modes >= margin
    if id.energyPD
        nonl = @(x) monoNonlcon(x, Uhx, Vbx, cI, del);
    else                                                % only Cbar modes >= margin
        Alin = [Alin; -Vbx];  blin = [blin; -mrg*ones(3*nK,1)];
        nonl = [];
    end
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
            'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
            'MaxIterations', 5000, 'MaxFunctionEvaluations', 1e6, ...
            'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-15);
        [x, ~, flag] = fmincon(obj, x0, Alin, blin, [], [], lb, ub, nonl, opts);
        method = sprintf('fmincon (sqp), exitflag %d', flag);
    elseif exist('sqp', 'file') == 2 || exist('sqp', 'builtin') == 5
        if isempty(nonl)
            hf = {@(x) blin - Alin*x, @(x) -Alin};
        else
            hf = {@(x) [blin - Alin*x; -firstOut(nonl, x)], @(x) [-Alin; -thirdOutT(nonl, x)]};
        end
        [x, ~, info] = sqp(x0, {@(x) firstOut(obj, x), @(x) secondOut(obj, x)}, [], hf, lb, ub, 2000, 1e-12);
        method = sprintf('sqp (Octave), info %d', info);
    else
        error('solveMonotone needs fmincon (Optimization Toolbox) or Octave''s sqp.');
    end

    c    = sc(:) .* x(:);
    coef = reshape(c, nb, np);                         % coef(k, j): parameter j, exponent e(k)
    coef(coef < 1e-10*max(abs(coef(:)))) = 0;          % round-off at the bound c >= 0
    P    = Bn * coef;                                  % nK x np
    res  = Pb*coef(:) - bb;
    rms  = nan(nK, 3);
    for iN = 1:nK
        for q = 1:3
            mq = use & chan == q & nIdx == iN;
            if any(mq), rms(iN,q) = sqrt(mean((res(mq)/sqrt(wq(q))).^2)); end
        end
    end
    mineig = zeros(nK, 3);
    for iN = 1:nK, mineig(iN,:) = minEigs(P(iN,:), C); end
    constantInN = pNames(all(coef(e > 0, :) == 0, 1));
    M = struct('params', P, 'coef', coef, 'exponents', e, 'Nlist', Nk.', 'rmsPerChannel', rms, ...
               'minEig', mineig, 'method', method, 'constantInN', {constantInN}, ...
               'objective', 0.5*norm(res(use))^2);
    M.eval = @(Nq) (Nq(:) .^ (-e)) * coef;
end

function [J, g] = monoObjective(x, P, b)
    r = P*x(:) - b;
    J = 0.5*(r.'*r);
    if nargout > 1, g = P.'*r; end
end

function [c, ceq, gc, gceq] = monoNonlcon(x, U, V, cI, del)
    % chat_i^2 / ((1-del) c_i) - cbar_i <= 0   (convex), one per mode and N
    uc = U*x(:);  vc = V*x(:);  k = 1 ./ ((1 - del) * cI);
    c   = k .* uc.^2 - vc;
    ceq = [];
    if nargout > 2
        gc = (2*(k .* uc) .* U - V).';                 % n x ncon
        gceq = [];
    end
end

function y = firstOut(f, x)
    y = f(x);
end

function y = secondOut(f, x)
    [~, y] = f(x);
end

function y = thirdOutT(f, x)
    [~, ~, y] = f(x);  y = y.';
end

function printMonotone(M, C, lamC, pN)
    Nk = M.Nlist(:);  P = M.params;
    fprintf('\n[identifyCbarChatG] monotone in N, all N together (%s)\n', M.method);
    fprintf('  C_hom: C11 = %.6g, C12 = %.6g, C22 = %.6g, C33 = %.6g\n', C(1,1), C(1,2), C(2,2), C(3,3));
    fprintf('  %4s', 'N');  fprintf(' %11s', pN{:});
    fprintf(' | %8s %8s %8s | %7s %7s %7s\n', 'eigCb', 'eigCh', 'eigTot', 'rmsE', 'rmsS', 'rmsD');
    for iN = 1:numel(Nk)
        fprintf('  %4d', Nk(iN));  fprintf(' %11.5g', P(iN,:));
        fprintf(' | %8.3g %8.3g %8.3g | %7.4f %7.4f %7.4f\n', M.minEig(iN,:)/lamC, M.rmsPerChannel(iN,:));
    end
    fprintf('  (eig: smallest eigenvalue of Cbar, Chat, [C Chat; Chat Cbar] / max eig(C) -- all >= 0)\n');
    fprintf('  normalized by C_hom:   N   Cb11/C11  Cb12/C12  Cb33/C33 | Ch11/C11  Ch12/C12  Ch33/C33\n');
    for iN = 1:numel(Nk)
        fprintf('                      %4d %9.4f %9.4f %9.4f | %9.4f %9.4f %9.4f\n', Nk(iN), ...
            P(iN,1)/C(1,1), P(iN,2)/C(1,2), P(iN,3)/C(3,3), P(iN,4)/C(1,1), P(iN,5)/C(1,2), P(iN,6)/C(3,3));
    end
    fprintf('  p(N) = sum_k c_k N^-e_k, c_k >= 0, e = %s:\n', mat2str(M.exponents));
    for j = 1:numel(pN)
        parts = arrayfun(@(k) termText(M.coef(k,j), M.exponents(k)), 1:numel(M.exponents), 'UniformOutput', false);
        fprintf('    %-5s(N) = %s\n', pN{j}, strtrim(strjoin(parts, ' ')));
    end
    if ~isempty(M.constantInN)
        fprintf(['  constant in N (the data would let them grow with N; the monotone-decrease\n' ...
                 '  constraint keeps them constant): %s\n'], strjoin(M.constantInN, ', '));
    end
end

function plotMonotone(R, C, lamC)
    M = R.monotone;  N = M.Nlist(:);  P = M.params;  Nf = linspace(min(N), max(N), 300).';
    Pf = M.eval(Nf);
    nrm = [C(1,1) C(1,2) C(3,3) C(1,1) C(1,2) C(3,3)];
    free = isfield(R, 'perN');
    figure('Color','w', 'Name', 'Params_vsN_monotone', 'Position', [40 60 1500 850]);
    grp = {1:3, 4:6};  ttl = {'(a) Cbar / C_{hom}', '(b) Chat / C_{hom}'};
    lb  = {{'Cbar_{11}/C_{11}', 'Cbar_{12}/C_{12}', 'Cbar_{33}/C_{33}'}, ...
           {'Chat_{11}/C_{11}', 'Chat_{12}/C_{12}', 'Chat_{33}/C_{33}'}};
    mk = {'o','s','^'};  cl = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
    for g = 1:2
        ax = subplot(2,2,g);  hold(ax,'on');  grid(ax,'on');  hh = [];
        for q = 1:3
            j = grp{g}(q);
            hh(end+1) = plot(ax, N, P(:,j)/nrm(j), mk{q}, 'Color', cl(q,:), 'MarkerFaceColor', cl(q,:)); %#ok<AGROW>
            plot(ax, Nf, Pf(:,j)/nrm(j), '-', 'Color', cl(q,:), 'LineWidth', 1.3);
            if free, plot(ax, N, R.perN.params(:,j)/nrm(j), mk{q}, 'Color', cl(q,:)); end
        end
        legend(ax, hh, lb{g}, 'Location', 'best');
        set(ax, 'XTick', N);  xlabel(ax, 'N');  title(ax, ttl{g});
    end
    ax = subplot(2,2,3);  hold(ax,'on');  grid(ax,'on');  hh = [];  lg = {};
    for j = 7:numel(R.paramNames)
        hh(end+1) = plot(ax, N, P(:,j), mk{j-6}, 'Color', cl(j-6,:), 'MarkerFaceColor', cl(j-6,:)); %#ok<AGROW>
        plot(ax, Nf, Pf(:,j), '-', 'Color', cl(j-6,:), 'LineWidth', 1.3);
        if free, plot(ax, N, R.perN.params(:,j), mk{j-6}, 'Color', cl(j-6,:)); end
        lg{end+1} = R.paramNames{j}; %#ok<AGROW>
    end
    legend(ax, hh, lg, 'Location', 'best');  set(ax, 'XTick', N);  xlabel(ax, 'N');
    title(ax, '(c) G');
    ax = subplot(2,2,4);  hold(ax,'on');  grid(ax,'on');
    plot(ax, N, M.minEig/lamC, '-o');  plot(ax, [min(N) max(N)], [0 0], 'k:');
    legend(ax, {'min eig Cbar', 'min eig Chat', 'min eig [C Chat; Chat Cbar]'}, 'Location', 'best');
    set(ax, 'XTick', N);  xlabel(ax, 'N');  ylabel(ax, '/ max eig(C_{hom})');
    title(ax, '(d) positive definiteness');
    if free
        sgtitle('Monotone Cbar(N), Chat(N), G(N): filled = identified, line = p(N), open = each N on its own');
    else
        sgtitle('Monotone Cbar(N), Chat(N), G(N)');
    end
end

function [Gc, h, names] = pdConstraints(nPar, id, lamC, pNames)
    % rows of Gc*theta >= h (theta = [Cb11 Cb12 Cb33 Ch11 Ch12 Ch33 G...])
    m = id.pdMargin * lamC;
    Gc = zeros(0, nPar);  h = zeros(0,1);  names = {};
    if id.pdCbar
        Gc(end+1,1:3) = [1 -1 0];  Gc(end+1,1:3) = [1 1 0];  Gc(end+1,1:3) = [0 0 1];
        h = [h; m; m; m];  names = [names, {'Cb11-Cb12', 'Cb11+Cb12', 'Cb33'}];
    end
    if id.pdChat
        Gc(end+1,4:6) = [1 -1 0];  Gc(end+1,4:6) = [1 1 0];  Gc(end+1,4:6) = [0 0 1];
        h = [h; m; m; m];  names = [names, {'Ch11-Ch12', 'Ch11+Ch12', 'Ch33'}];
    end
    for j = 7:nPar                                        % G >= 0
        Gc(end+1, j) = 1;  h = [h; 0];  names{end+1} = pNames{j};   %#ok<AGROW>
    end
end

function [x, active, method] = qpActiveSet(P, b, G, h, names)
    % min ||P x - b||^2  s.t.  G x >= h, solved exactly: the optimum is the
    % equality-constrained least-squares solution of its own active set, so
    % the best feasible candidate over all active sets is the global optimum
    H = P.'*P;  f = P.'*b;  n = size(P,2);  nc = size(G,1);
    if nc > 0                                   % balance the KKT system: unit constraint rows
        rn = sqrt(sum(G.^2, 2));  G = G ./ rn;  h = h ./ rn;
    end
    tol = 1e-9 * max(1, max(abs(h)));
    x = pinv(H)*f;
    if isempty(G) || all(G*x >= h - tol)
        active = {};  method = 'unconstrained least squares';  return;
    end
    best = inf;  xb = x;  ab = [];
    for mask = 1:(2^nc - 1)
        Sa = find(bitget(mask, 1:nc));
        Ga = G(Sa,:);
        K  = [H, Ga.'; Ga, zeros(numel(Sa))];
        z  = pinv(K) * [f; h(Sa)];
        xc = z(1:n);
        if all(G*xc >= h - tol)
            J = norm(P*xc - b)^2;
            if J < best, best = J;  xb = xc;  ab = Sa; end
        end
    end
    if isinf(best)
        error('identifyCbarChatG: no feasible point for the positive-definiteness constraints.');
    end
    x = xb;
    active = names(ab);
    method = 'least squares, PD bounds active';
end

function [theta, method, active] = energyPDsolve(theta0, P, b, Gc, h, C, id)
    % additionally keep [C Chat; Chat Cbar] positive definite (non-linear)
    m = id.pdMargin * max(eig(C));
    obj = @(t) sum((P*t(:) - b).^2);
    nonl = @(t) m - min(eig(blockMatrix(t, C)));          % <= 0
    active = {};
    if -nonl(theta0) >= 0
        theta = theta0;  method = 'least squares (energy block already positive definite)';  return;
    end
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIterations', 3000);
        theta = fmincon(@(t) obj(t), theta0(:), -Gc, -h, [], [], [], [], ...
                        @(t) pdNonlcon(t, C, m), opts).';
        method = 'fmincon, energy block PD';
    elseif exist('sqp', 'file') == 2 || exist('sqp', 'builtin') == 5
        theta = sqp(theta0(:), obj, [], @(t) [Gc*t(:) - h; -nonl(t)], [], [], 500).';
        method = 'sqp (Octave), energy block PD';
    else
        theta = theta0;  method = 'energy block PD not enforced (no solver)';
    end
    active = {'energy block PD'};
end

function [c, ceq] = pdNonlcon(t, C, m)
    c = m - min(eig(blockMatrix(t, C)));
    ceq = [];
end

function B = blockMatrix(t, C)
    Cb = [t(1) t(2) 0; t(2) t(1) 0; 0 0 t(3)];
    Ch = [t(4) t(5) 0; t(5) t(4) 0; 0 0 t(6)];
    B  = [C Ch; Ch Cb];
end

function e = minEigs(t, C)
    Cb = [t(1) t(2) 0; t(2) t(1) 0; 0 0 t(3)];
    Ch = [t(4) t(5) 0; t(5) t(4) 0; 0 0 t(6)];
    e  = [min(eig(Cb)), min(eig(Ch)), min(eig(blockMatrix(t, C)))];
end

% ============================================================
% ================= trend, report, output ====================
% ============================================================

function T = fitTrends(N, P, names)
    N = N(:);  nK = numel(N);
    if nK >= 5, e = [0 1 2]; elseif nK >= 3, e = [0 1]; else, e = 0; end
    B = N .^ (-e);
    T.exponents = e;  T.coef = B \ P;  T.rms = sqrt(mean((B*T.coef - P).^2, 1));
    T.text = cell(1, numel(names));
    for j = 1:numel(names)
        parts = arrayfun(@(k) termText(T.coef(k,j), e(k)), 1:numel(e), 'UniformOutput', false);
        T.text{j} = strtrim(strjoin(parts, ' '));
    end
    T.eval = @(Nq, j) (Nq(:) .^ (-e)) * T.coef(:,j);
end

function s = termText(c, p)
    if p == 0, s = sprintf('%.6g', c);
    elseif p == 1, s = sprintf('%+.6g/N', c);
    else, s = sprintf('%+.6g/N^%d', c, p);
    end
end

function printPerN(R, C, lamC)
    Nk = R.Nlist(:);  P = R.perN.params;  pN = R.paramNames;
    fprintf('\n[identifyCbarChatG] parameters for each N (%d load cases per N)\n', numel(R.cases));
    fprintf('  C_hom: C11 = %.6g, C12 = %.6g, C22 = %.6g, C33 = %.6g\n', C(1,1), C(1,2), C(2,2), C(3,3));
    fprintf('  %4s', 'N');  fprintf(' %11s', pN{:});
    fprintf(' | %8s %8s %8s | %7s %7s %7s | %s\n', 'eigCb', 'eigCh', 'eigTot', 'rmsE', 'rmsS', 'rmsD', 'active PD bounds');
    for iN = 1:numel(Nk)
        fprintf('  %4d', Nk(iN));  fprintf(' %11.5g', P(iN,:));
        fprintf(' | %8.3g %8.3g %8.3g | %7.4f %7.4f %7.4f | %s\n', R.perN.minEig(iN,:)/lamC, ...
            R.perN.rmsPerChannel(iN,:), strjoin(R.perN.activeConstraints{iN}, ', '));
    end
    fprintf('  (eig columns: smallest eigenvalue of Cbar, Chat, [C Chat; Chat Cbar], divided by max eig(C))\n');
    fprintf('  normalized by C_hom:   N   Cb11/C11  Cb12/C12  Cb33/C33 | Ch11/C11  Ch12/C12  Ch33/C33\n');
    for iN = 1:numel(Nk)
        fprintf('                      %4d %9.4f %9.4f %9.4f | %9.4f %9.4f %9.4f\n', Nk(iN), ...
            P(iN,1)/C(1,1), P(iN,2)/C(1,2), P(iN,3)/C(3,3), P(iN,4)/C(1,1), P(iN,5)/C(1,2), P(iN,6)/C(3,3));
    end
    fprintf('  informative load cases per channel at N = %d (energy/stress/double): %s\n', ...
        Nk(end), mat2str(R.perN.nInformative(end,:)));
    fprintf('  N-dependence of the free per-N values, p ~ N^s (s = 0: material constant):\n    ');
    for j = 1:numel(pN)
        v = P(:,j);
        if all(v > 0), sl = [log(Nk(:)), ones(numel(Nk),1)] \ log(v); fprintf('%s: s=%+.2f   ', pN{j}, sl(1));
        else, fprintf('%s: (sign change)   ', pN{j}); end
    end
    fprintf('\n');
    fprintf('  trend  p(N) = c0 + c1/N%s:\n', ternary(numel(R.trend.exponents) > 2, ' + c2/N^2', ''));
    for j = 1:numel(pN)
        fprintf('    %-5s = %s   (limit N->inf %.6g)\n', pN{j}, R.trend.text{j}, R.trend.coef(1,j));
    end
end

function writeCSV(fname, R)
    fid = fopen(fname, 'wt');  if fid < 0, return; end
    if isfield(R, 'anchored')
        M = R.anchored;
        fprintf(fid, 'N');  fprintf(fid, ',%s', R.paramNames{:});
        fprintf(fid, ',minEig_Cbar,minEig_Chat,minEig_block,rms_energy,rms_stress,rms_double\n');
        for iN = 1:numel(M.Nlist)
            fprintf(fid, '%d', M.Nlist(iN));
            fprintf(fid, ',%.10g', M.params(iN,:), M.minEig(iN,:), M.rmsPerChannel(iN,:));
            fprintf(fid, '\n');
        end
        fprintf(fid, 'limit N->inf');  fprintf(fid, ',%.10g', M.anchor(1:6), M.G);  fprintf(fid, '\n');
        for k = 1:numel(M.exponents)
            fprintf(fid, 'deviation N^-%g', M.exponents(k));  fprintf(fid, ',%.10g', M.Dev(k,:));
            fprintf(fid, '\n');
        end
        fprintf(fid, 'm,%.10g\n', M.m);
        fclose(fid);  return;
    end
    if isfield(R, 'monotone')
        M = R.monotone;
        fprintf(fid, 'N');  fprintf(fid, ',%s', R.paramNames{:});
        fprintf(fid, ',minEig_Cbar,minEig_Chat,minEig_block,rms_energy,rms_stress,rms_double\n');
        for iN = 1:numel(M.Nlist)
            fprintf(fid, '%d', M.Nlist(iN));
            fprintf(fid, ',%.10g', M.params(iN,:), M.minEig(iN,:), M.rmsPerChannel(iN,:));
            fprintf(fid, '\n');
        end
        for k = 1:numel(M.exponents)
            fprintf(fid, 'coef N^-%g', M.exponents(k));  fprintf(fid, ',%.10g', M.coef(k,:));
            fprintf(fid, '\n');
        end
        fclose(fid);  return;
    end
    fprintf(fid, 'N');  fprintf(fid, ',%s', R.paramNames{:});
    fprintf(fid, ',minEig_Cbar,minEig_Chat,minEig_block,rms_energy,rms_stress,rms_double\n');
    if isfield(R, 'perN')
        for iN = 1:numel(R.Nlist)
            fprintf(fid, '%d', R.Nlist(iN));
            fprintf(fid, ',%.10g', R.perN.params(iN,:), R.perN.minEig(iN,:), R.perN.rmsPerChannel(iN,:));
            fprintf(fid, '\n');
        end
    end
    if isfield(R, 'global')
        fprintf(fid, 'all');
        fprintf(fid, ',%.10g', R.global.params, R.global.minEig, R.global.rmsPerChannel);
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function plotParams(R, C, lamC)
    N = R.Nlist(:);  P = R.perN.params;  Nf = linspace(min(N), max(N), 200).';
    nrm = [C(1,1) C(1,2) C(3,3) C(1,1) C(1,2) C(3,3)];
    figure('Color','w', 'Name', 'Params_vsN', 'Position', [40 60 1500 850]);
    grp = {1:3, 4:6};  ttl = {'(a) Cbar / C_{hom}', '(b) Chat / C_{hom}'};
    lb  = {{'Cbar_{11}/C_{11}', 'Cbar_{12}/C_{12}', 'Cbar_{33}/C_{33}'}, ...
           {'Chat_{11}/C_{11}', 'Chat_{12}/C_{12}', 'Chat_{33}/C_{33}'}};
    mk = {'o','s','^'};  cl = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
    for g = 1:2
        ax = subplot(2,2,g);  hold(ax,'on');  grid(ax,'on');  hh = [];
        for q = 1:3
            j = grp{g}(q);
            hh(end+1) = plot(ax, N, P(:,j)/nrm(j), mk{q}, 'Color', cl(q,:), 'MarkerFaceColor', cl(q,:)); %#ok<AGROW>
            plot(ax, Nf, R.trend.eval(Nf, j)/nrm(j), '-', 'Color', cl(q,:));
            if isfield(R, 'global'), plot(ax, [min(N) max(N)], R.global.params(j)/nrm(j)*[1 1], '--', 'Color', cl(q,:)); end
        end
        legend(ax, hh, lb{g}, 'Location', 'best', 'Interpreter', 'tex');
        set(ax, 'XTick', N);  xlabel(ax, 'N');  title(ax, ttl{g}, 'Interpreter', 'tex');
    end
    ax = subplot(2,2,3);  hold(ax,'on');  grid(ax,'on');  hh = [];  lg = {};
    for j = 7:numel(R.paramNames)
        hh(end+1) = plot(ax, N, P(:,j), mk{j-6}, 'Color', cl(j-6,:), 'MarkerFaceColor', cl(j-6,:)); %#ok<AGROW>
        plot(ax, Nf, R.trend.eval(Nf, j), '-', 'Color', cl(j-6,:));
        lg{end+1} = R.paramNames{j}; %#ok<AGROW>
    end
    legend(ax, hh, lg, 'Location', 'best');  set(ax, 'XTick', N);  xlabel(ax, 'N');
    title(ax, '(c) G  (trend line; dashed = one set for all N)');
    ax = subplot(2,2,4);  hold(ax,'on');  grid(ax,'on');
    plot(ax, N, R.perN.minEig/lamC, '-o');
    plot(ax, [min(N) max(N)], [0 0], 'k:');
    legend(ax, {'min eig Cbar', 'min eig Chat', 'min eig [C Chat; Chat Cbar]'}, 'Location', 'best');
    set(ax, 'XTick', N);  xlabel(ax, 'N');  ylabel(ax, '/ max eig(C_{hom})');
    title(ax, '(d) positive definiteness');
    sgtitle('Identified Cbar, Chat, G as functions of N');
end

function plotFit(S, Pn, cases, id, cfg, C, figName, ttl)
    Nk = S.Nlist(:);  nC = numel(cases);  cm = 0.85*hsv(nC);
    idu = id;  idu.wE = 1;  idu.wS = 1;  idu.wD = 1;       % unweighted residuals for plotting
    figure('Color','w', 'Name', figName, 'Position', [30 30 1650 900]);
    axE = subplot(2,3,1); hold(axE,'on'); grid(axE,'on');
    axS = subplot(2,3,2); hold(axS,'on'); grid(axS,'on');
    axM = subplot(2,3,3); hold(axM,'on'); grid(axM,'on');
    axSc = subplot(2,3,4); hold(axSc,'on'); grid(axSc,'on');
    axMc = subplot(2,3,5); hold(axMc,'on'); grid(axMc,'on');
    for t = 1:nC
        ic = cases(t);  eE = nan(numel(Nk),1);  eS = eE;  eM = eE;
        for iN = 1:numel(Nk)
            [Pr, br, ch] = caseRows(S, ic, iN, idu, cfg, C);
            res = Pr*Pn(iN,:).' - br;                          % scaled residuals
            eE(iN) = res(ch == 1);
            if any(ch == 2)
                eS(iN) = norm(res(ch == 2));
                fe = br(ch == 2);  mo = Pr(ch == 2,:)*Pn(iN,:).';
                plot(axSc, fe, mo, '.', 'Color', cm(t,:), 'MarkerSize', 8);
            end
            if any(ch == 3)
                eM(iN) = norm(res(ch == 3));
                fe = br(ch == 3);  mo = Pr(ch == 3,:)*Pn(iN,:).';
                plot(axMc, fe, mo, '.', 'Color', cm(t,:), 'MarkerSize', 8);
            end
        end
        plot(axE, Nk, eE, '-o', 'Color', cm(t,:), 'MarkerSize', 4);
        plot(axS, Nk, eS, '-o', 'Color', cm(t,:), 'MarkerSize', 4);
        plot(axM, Nk, eM, '-o', 'Color', cm(t,:), 'MarkerSize', 4);
    end
    title(axE, '(a) energy: (W_{model} - 2U_{het}) / 2U_{het}');
    title(axS, '(b) stress: |model - FE| / |FE|');
    title(axM, '(c) double stress: |model - FE| / |FE|');
    for ax = [axE axS axM], set(ax, 'XTick', Nk); xlabel(ax, 'N'); end
    for ax = [axSc axMc]
        lim = [min([get(ax,'XLim') get(ax,'YLim')]), max([get(ax,'XLim') get(ax,'YLim')])];
        plot(ax, lim, lim, 'k--');  xlabel(ax, 'FE (scaled)');  ylabel(ax, 'model (scaled)');  axis(ax, 'square');
    end
    title(axSc, '(d) \int\sigma components (scaled by |FE|)');
    title(axMc, '(e) \int\sigma_{ij}x_k components (scaled by |FE|)');
    sgtitle(ttl);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
