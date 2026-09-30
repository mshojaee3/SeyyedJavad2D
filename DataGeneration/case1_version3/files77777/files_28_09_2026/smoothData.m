function S = smoothData(D, cfg)
% SMOOTHDATA  Remove the numerical fluctuation over N from the FE data.
%
%   S = smoothData(D, cfg)      D from dataGeneration.m
%
% Every quantity is fitted ON ITS OWN over N >= cfg.smooth.NminFit (no
% relation between the integrals is imposed), by y(N) ~ sum_k c_k N^-p_k
% with fewer terms than values of N:
%
%   quantity                              model                       behaviour
%   A5, A7, A8        gamma squares       c_k >= 0, p in decayExp     decreasing -> 0
%   A6, A9..A13       gamma cross terms   one sign, p in decayExp     |.| decreasing -> 0
%   I19..I24          int gamma x_k       one sign, p in decayExp     |.| decreasing -> 0
%   A14..A19          kappa squares       ell^2*A fitted: c_k >= 0,   ell^2 A decreasing -> 0
%                                         p in decayExp               (bulk ~N^-2, boundary layer ~N^-1;
%                                                                     A itself may grow with N)
%   2U_het            2U_hom + one-sign decay of the excess          excess -> 0
%   A1..A4, I1..I18 (except I19..I24), Sfe, Mfe
%                     free signs, p in freeExp (0 1 2)                smooth
%
% Weighted relative least squares, w = 1/max(|y|, 1e-3 max|y|).
% All quantities are fitted as DENSITIES: divided by the region volume V(N)
% (moments I13..I24 and int sigma x additionally by the rms moment arm Lr(N)),
% then multiplied back. For the full domain V and Lr are constant; for the
% interior region they change with N, and this geometric N-dependence must not
% be mistaken for a trend of the data.

if ~isfield(D, 'A') || ~isfield(D, 'I') || ~isfield(D, 'VolN')
    error(['smoothData: DataGeneration.mat is from an older version. ' ...
           'Run dataGeneration again (cfg.runDataGeneration = true).']);
end
region = 'full';  if isfield(cfg, 'region'), region = cfg.region; end
D = selectRegion(D, region);
sm   = cfg.smooth;
Nl   = D.Nlist(:);
keep = Nl >= sm.NminFit & D.valid(:);
if nnz(keep) < 3
    error('smoothData: region ''%s'' has only %d values of N >= NminFit (need 3).', region, nnz(keep));
end
Nk = Nl(keep);  nK = numel(Nk);  nC = size(D.twoU_het, 2);
Pd = sm.decayExponents;  P4 = sm.T4Exponents;  Pf = sm.freeExponents;

modeA = repmat({'free'}, 1, 19);   expA = repmat({Pf}, 1, 19);
modeA([5 7 8]) = {'positive'};     expA([5 7 8]) = {Pd};
modeA([6 9:13]) = {'oneSign'};     expA([6 9:13]) = {Pd};
modeA(14:19) = {'positive'};       expA(14:19) = {Pd};   % applied to ell^2*A
modeI = repmat({'free'}, 1, 24);   expI = repmat({Pf}, 1, 24);
modeI(19:24) = {'oneSign'};        expI(19:24) = {Pd};

S = struct();
S.A   = zeros(nK, nC, 19);  S.I = zeros(nK, nC, 24);
S.Sfe = zeros(nK, nC, 3);   S.Mfe = zeros(nK, nC, 6);
S.twoU_het = zeros(nK, nC); S.twoU_hom = zeros(nK, nC);
rmsA = zeros(nC, 19);  rmsI = zeros(nC, 24);  rmsS = zeros(nC, 3);  rmsM = zeros(nC, 6);  rmsU = zeros(nC,1);

V  = D.VolN(keep);  V = V(:);                 % region volume per N
l2 = D.ell2(keep);  l2 = l2(:);               % (L/N)^2
VL = V .* D.LrN(keep);  VL = VL(:);           % volume x rms moment arm
sI = [repmat(V,1,12), repmat(VL,1,12)];       % I1..I12 volume, I13..I24 moments
for ic = 1:nC
    for j = 1:19
        if j >= 14, sc = V ./ l2; else, sc = V; end      % kappa terms: fit ell^2*A/V
        [y, ~, rmsA(ic,j)] = fitN(Nk, D.A(keep,ic,j)./sc, expA{j}, modeA{j});
        S.A(:,ic,j) = y .* sc;
    end
    for j = 1:24
        [y, ~, rmsI(ic,j)] = fitN(Nk, D.I(keep,ic,j)./sI(:,j), expI{j}, modeI{j});
        S.I(:,ic,j) = y .* sI(:,j);
    end
    for j = 1:3
        [y, ~, rmsS(ic,j)] = fitN(Nk, D.Sfe(keep,ic,j)./V, Pf, 'free');   S.Sfe(:,ic,j) = y .* V;
    end
    for j = 1:6
        [y, ~, rmsM(ic,j)] = fitN(Nk, D.Mfe(keep,ic,j)./VL, Pf, 'free');  S.Mfe(:,ic,j) = y .* VL;
    end
    eh = D.twoU_hom(keep,ic) ./ V;             % Stage-1 energy density (exact, per N)
    [dU, ~, rmsU(ic)] = fitN(Nk, D.twoU_het(keep,ic)./V - eh, Pd, 'oneSign');
    S.twoU_hom(:,ic) = D.twoU_hom(keep,ic);
    S.twoU_het(:,ic) = (eh + dU) .* V;
end

S.Nlist = Nk.';  S.ell2 = D.ell2(keep);  S.NminFit = sm.NminFit;  S.keep = keep.';
S.caseNames = D.caseNames;  S.caseLabels = D.caseLabels;  S.rawVecs = D.rawVecs;
S.C_hom = D.C_hom;  S.VolN = D.VolN(keep);  S.LrN = D.LrN(keep);  S.region = region;  S.raw = D;
S.fitRMS = struct('A', rmsA, 'I', rmsI, 'Sfe', rmsS, 'Mfe', rmsM, 'twoU', rmsU);

% ---- report: worst relative RMS per group of quantities ----
grp = {'psi.psi   A1-A4',   rmsA(:,1:4);   'gam.gam   A5-A8',   rmsA(:,5:8);
       'psi.gam   A9-A13',  rmsA(:,9:13);  'kap.kap   A14-A19', rmsA(:,14:19);
       'int psi   I1-I3',   rmsI(:,1:3);   'int gam   I4-I6',   rmsI(:,4:6);
       'int kap   I7-I12',  rmsI(:,7:12);  'psi x     I13-I18', rmsI(:,13:18);
       'gam x     I19-I24', rmsI(:,19:24); 'int sigma (FE)',     rmsS;
       'int sigma x (FE)',  rmsM;          '2U_het',             rmsU};
fprintf('[smoothData] region ''%s'', N = %s, %d cases, every quantity fitted on its own\n', region, mat2str(Nk.'), nC);
fprintf('  %-20s %12s %12s\n', 'group', 'median RMS', 'max RMS');
for g = 1:size(grp,1)
    v = grp{g,2}(:);
    fprintf('  %-20s %12.3g %12.3g\n', grp{g,1}, median(v), max(v));
end
fprintf('  (relative RMS of the fits; quantities that are ~0 for a load case by symmetry show large relative values)\n');

save(fullfile(cfg.outDir, 'SmoothedData.mat'), '-struct', 'S');
if cfg.plots, plotSmoothing(D, S, keep); end
end


% ============================================================

function D = selectRegion(D, region)
    % map the chosen region onto the standard field names
    switch lower(region)
        case 'full'
            D.valid = true(numel(D.Nlist), 1);
        case 'interior'
            if ~isfield(D, 'int') || ~any(D.int.valid)
                error('smoothData: no interior data (set cfg.interiorBand > 0 and rerun dataGeneration).');
            end
            for f = {'A','AM','I','IM','Sfe','Mfe','twoU_het','twoU_het_std','twoU_hom','VolN','LrN'}
                if isfield(D.int, f{1}), D.(f{1}) = D.int.(f{1}); end
            end
            D.valid = D.int.valid(:);
        otherwise
            error('smoothData: unknown region ''%s''.', region);
    end
end

function [ys, c, rmsRel] = fitN(N, y, P, mode)
    % y(N) ~ sum_k c_k N^-P(k)
    %   'positive'       : c >= 0
    %   'offsetPositive' : constant free, decay terms c >= 0
    %   'oneSign'        : c >= 0 or c <= 0, whichever fits better
    %   'free'           : unconstrained
    N = N(:);  y = y(:);  P = P(:).';
    if strcmp(mode, 'offsetPositive'), P = P(P > 0); nMax = numel(N) - 2; else, nMax = numel(N) - 1; end
    nMax = max(nMax, 1);
    if numel(P) > nMax                           % fewer terms than points: best subset of exponents
        sub = nchoosek(1:numel(P), nMax);  best = inf;
        for q = 1:size(sub,1)
            [ysq, cq, rq] = fitN(N, y, P(sub(q,:)), mode);
            if rq < best, best = rq;  ys = ysq;  c = cq;  rmsRel = rq; end
        end
        return;
    end
    A = N .^ (-P);
    sc = max(abs(y));  if sc == 0, ys = zeros(size(y)); c = zeros(numel(P),1); rmsRel = 0; return; end
    w = 1 ./ max(abs(y), 1e-3*sc);               % relative weighting
    Aw = A .* w;  yw = y .* w;
    switch mode
        case 'positive'
            c = nnlsScaled(Aw, yw);
        case 'offsetPositive'
            cc = nnlsScaled([w, -w, Aw], yw);
            c  = [cc(1) - cc(2); cc(3:end)];
            A  = [ones(size(N)), A];
        case 'oneSign'
            cp = nnlsScaled(Aw,  yw);  cm = nnlsScaled(Aw, -yw);
            if norm(Aw*cp - yw) <= norm(Aw*cm + yw), c = cp; else, c = -cm; end
        case 'free'
            c = Aw \ yw;
        otherwise
            error('fitN: unknown mode %s', mode);
    end
    ys = A * c;
    rmsRel = sqrt(mean(((ys - y) ./ max(abs(y), 1e-3*sc)).^2));
end

function c = nnlsScaled(A, y)
    % non-negative least squares on column- and row-balanced data
    % (avoids the LSQNONNEG iteration limit on badly scaled columns)
    sy = max(abs(y));  if sy == 0, c = zeros(size(A,2),1); return; end
    d  = max(abs(A), [], 1);  d(d == 0) = 1;
    c  = lsqnonneg(A ./ d, y / sy);
    c  = c(:) ./ d(:) * sy;
end

function plotSmoothing(D, S, keep)
    Nl = D.Nlist(:);  Nk = S.Nlist(:);  nC = size(S.A, 2);  cm = 0.85*hsv(nC);
    C  = S.C_hom;
    % gamma-energy-like and coupling-like combinations (for display only)
    gE = @(X) X(:,:,5)*C(1,1) + 2*X(:,:,6)*C(1,2) + X(:,:,7)*C(2,2) + 4*X(:,:,8)*C(3,3);
    pE = @(X) X(:,:,1)*C(1,1) + 2*X(:,:,2)*C(1,2) + X(:,:,3)*C(2,2) + 4*X(:,:,4)*C(3,3);
    cE = @(X) (X(:,:,9)+X(:,:,12))*C(1,1) + (X(:,:,10)+X(:,:,11))*C(1,2) + 4*X(:,:,13)*C(3,3);
    kS = @(X) sum(X(:,:,14:19), 3);
    items = {'\int\gamma^TC\gamma / 2U_{hom}', gE;  '\int\psi^TC\psi / 2U_{hom}', pE;
             '\int\psi^TC\gamma / 2U_{hom}', cE;    '\int\Sigma\kappa^2 (A14..A19)', kS};
    figure('Color','w', 'Name', 'Smoothing', 'Position', [30 30 1750 950]);
    for q = 1:6
        ax = subplot(2,3,q);  hold(ax,'on');  grid(ax,'on');
        for ic = 1:nC
            if q <= 4
                r = items{q,2}(D.A(:,ic,:));  s = items{q,2}(S.A(:,ic,:));
                if q <= 3, r = r ./ D.twoU_hom(:,ic);  s = s ./ S.twoU_hom(:,ic); end
            elseif q == 5
                r = D.twoU_het(:,ic)./D.twoU_hom(:,ic);  s = S.twoU_het(:,ic)./S.twoU_hom(:,ic);
            else
                r = sqrt(sum(D.Mfe(:,ic,:).^2, 3));  s = sqrt(sum(S.Mfe(:,ic,:).^2, 3));
            end
            plot(ax, Nl(keep), r(keep), 'o', 'Color', cm(ic,:), 'MarkerSize', 3);
            plot(ax, Nk, s, '-', 'Color', cm(ic,:), 'LineWidth', 1.2);
        end
        set(ax, 'XTick', Nk);  xlim(ax, [Nk(1)-0.2, Nk(end)+0.2]);  xlabel(ax, 'N');
        if q <= 4, title(ax, items{q,1});
        elseif q == 5, title(ax, '2U_{het} / 2U_{hom}');
        else, title(ax, '|\int\sigma_{ij}x_k| (FE)'); end
    end
    sgtitle('smoothData: raw (o) and smoothed (lines), every quantity fitted separately');
end
