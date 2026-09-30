function A = identifyAlphas(S, cfg)
% IDENTIFYALPHAS  Micromorphic parameters from energy, integrated stress and
%                 integrated double stress; one set of alphas for every N
%                 (and optionally one global set).
%
%   A = identifyAlphas(S, cfg)     S from smoothData.m
%
% Model (per unit volume, C = C_hom, ell = L/N):
%   W = am/2 psi'C psi + ah/2 gamma'C gamma + ac psi'C gamma + al/2 ell^2 kappa'A kappa
%   s_(sym)       = dW/dpsi   = am C psi + ac C gamma
%   sigma_(gamma) = dW/dgamma = ah C gamma + ac C psi   (conjugate to eps = sym grad u)
%   mu            = dW/dkappa = al ell^2 A kappa
%
% Equations per load case and N (T1..T5 used as independent data):
%   energy (1)        am T2 + ah T3 + 2 ac T5 + al ell^2 T4        = 2U_het
%   stress (3)        int sigma_(gamma) = ah S3 + ac S2             = Sfe = int sigma_FE
%     (+3 if cfg.id.stressModel = 'both':  int s_(sym) = am S2 + ac S3 = Sfe)
%   double stress (6) cfg.id.doubleStressModel
%     'local+mu' : int sigma_(gamma) x + int mu = ah M3 + ac M2 + al ell^2 S4 = Mfe
%     'mu'       : int mu = al ell^2 S4                                        = Mfe
%   with Sfe = int sigma dOmega and Mfe = int sigma_ij x_k dOmega of the
%   heterogeneous FE solution (mean over the inclusion shifts).
%
% Each equation group is divided by its size (2U_het, |Sfe|, |Mfe|) and
% weighted by cfg.id.wE / wS / wD. A stress or double-stress group is used
% only when it carries information: |Sfe| >= tolInfo * sigRef and
% |Mfe| >= tolInfo * sigRef * L/sqrt(12), sigRef = sqrt(2U_het*max eig(C)/V)*V
% (by symmetry int sigma = 0 for pure G loads and int sigma x = 0 for pure
% H loads -- such groups carry no information).
%
% cfg.id.mode: 'perN' (alphas for each N + trend c0 + c1/N + c2/N^2),
%              'global' (one set for all N), 'both' (default).
% Constraints: am fixed (cfg.id.fixAlphaM) or > 0, ah >= 0, al >= 0,
%              am*ah - ac^2 >= pdMargin*am*ah.

id = cfg.id;
if ~isfield(id, 'mode') || isempty(id.mode), id.mode = 'both'; end
nC = size(S.T1, 2);
cases = id.cases;  if isempty(cases), cases = 1:nC; end
Nk = S.Nlist(:);  nK = numel(Nk);
names = {'alpha_m', 'alpha_h', 'alpha_c', 'alpha_l'};
chName = {'energy', 'stress sigma_gamma', 'stress s_sym', 'double stress'};

A = struct('names', {names}, 'channelNames', {chName}, 'cases', cases, 'Nlist', Nk.', ...
           'weights', [id.wE id.wS id.wD], 'stressModel', id.stressModel, ...
           'doubleStressModel', id.doubleStressModel, 'mode', id.mode);

% ---------------- one set of alphas per N ----------------
if any(strcmpi(id.mode, {'perN', 'both'}))
    aN = nan(nK, 4);  rmsN = nan(nK, 4);  okN = false(nK, 1);  condN = nan(nK, 1);
    methN = cell(nK, 1);  estN = cell(nK, 1);  infoN = zeros(nK, 4);
    for iN = 1:nK
        r = solveBlock(S, cases, iN, id, cfg);
        aN(iN,:) = r.a;  rmsN(iN,:) = r.rmsCh;  okN(iN) = r.ok;  methN{iN} = r.method;
        condN(iN) = r.sv(1)/r.sv(end);  estN{iN} = r.est;  infoN(iN,:) = r.nInfo;
    end
    A.perN = struct('alphas', aN, 'rmsPerChannel', rmsN, 'feasible', okN, ...
                    'method', {methN}, 'cond', condN, 'channelEstimates', {estN}, 'nInformative', infoN);
    A.trend = fitTrends(Nk, aN, id);

    fprintf('\n[identifyAlphas] alphas for each N (%d load cases per N, wE=%g wS=%g wD=%g%s)\n', ...
        numel(cases), id.wE, id.wS, id.wD, ternary(id.fixAlphaM, sprintf(', alpha_m fixed = %g', id.alphaM), ''));
    fprintf('  %4s %11s %11s %11s %11s | %9s %9s | %7s %7s %7s %7s | %6s %s\n', 'N', names{:}, ...
        'am+ah-2ac', 'amah-ac^2', 'rmsE', 'rmsSg', 'rmsSs', 'rmsD', 'cond', 'PD');
    for iN = 1:nK
        a = aN(iN,:);
        fprintf('  %4d %11.5g %11.5g %11.5g %11.5g | %9.4g %9.4g | %7.4f %7.4f %7.4f %7.4f | %6.3g %s\n', ...
            Nk(iN), a, a(1)+a(2)-2*a(3), a(1)*a(2)-a(3)^2, rmsN(iN,:), condN(iN), ...
            ternary(okN(iN), 'ok', 'VIOLATED'));
    end
    onB = ~strcmp(methN, 'unconstrained least squares');
    fprintf('  N with the PD bound active (constrained solve): %s\n', mat2str(Nk(onB).'));
    fprintf('  informative load cases per channel (energy/stress/s_sym/double), N = %d: %s\n', ...
        Nk(end), mat2str(infoN(end,:)));
    fprintf('  channel-only estimates at N = %d [am ah ac al] (NaN = not seen by the channel):\n', Nk(end));
    for q = 1:4
        fprintf('    %-20s %s\n', chName{q}, mat2str(estN{end}(q,:), 5));
    end
    fprintf('  trend  alpha(N) = %s :\n', A.trend.formula);
    for j = 1:4
        if id.fixAlphaM && j == 1, continue; end
        fprintf('    %-8s = %s   (limit N->inf: %.5g, fit RMS %.2g)\n', names{j}, ...
            A.trend.text{j}, A.trend.coef(1,j), A.trend.rms(j));
    end
end

% ---------------- one set of alphas for all N ----------------
if any(strcmpi(id.mode, {'global', 'both'}))
    r = solveBlock(S, cases, 1:nK, id, cfg);
    A.global = struct('alphas', r.a, 'rmsPerChannel', r.rmsCh, 'feasible', r.ok, ...
                      'method', r.method, 'singularValues', r.sv, 'channelEstimates', r.est);
    fprintf('\n[identifyAlphas] one set for all N = %s: am=%.5g ah=%.5g ac=%.5g al=%.5g  (%s, PD %s)\n', ...
        mat2str(Nk.'), r.a, r.method, ternary(r.ok, 'ok', 'VIOLATED'));
    fprintf('  RMS misfit: energy %.4f | stress sigma_gamma %.4f | stress s_sym %.4f | double stress %.4f\n', r.rmsCh);
end

% convenience fields: per-N result if available (column over N)
if isfield(A, 'perN')
    A.alpha_m = A.perN.alphas(:,1);  A.alpha_h = A.perN.alphas(:,2);
    A.alpha_c = A.perN.alphas(:,3);  A.alpha_l = A.perN.alphas(:,4);
else
    A.alpha_m = A.global.alphas(1);  A.alpha_h = A.global.alphas(2);
    A.alpha_c = A.global.alphas(3);  A.alpha_l = A.global.alphas(4);
end

save(fullfile(cfg.outDir, 'IdentifiedAlphas.mat'), '-struct', 'A');
writeAlphaCSV(fullfile(cfg.outDir, 'IdentifiedAlphas_vsN.csv'), A);
if cfg.plots
    if isfield(A, 'perN')
        plotAlphasVsN(A, id);
        plotFit(S, A.perN.alphas, cases, 'IdentifiedAlphas_perN', ...
            sprintf('model with alpha(N) identified at each N, %d cases', numel(cases)), id, cfg);
    end
    if isfield(A, 'global')
        plotFit(S, repmat(A.global.alphas, nK, 1), cases, 'IdentifiedAlphas_global', ...
            sprintf('model with one alpha set for all N, %d cases', numel(cases)), id, cfg);
    end
end
end


% ============================================================
% ======================== solver ============================
% ============================================================

function r = solveBlock(S, cases, iNs, id, cfg)
    % weighted least squares over the given load cases and N indices
    V  = cfg.geom.Lx_tot * cfg.geom.Ly_tot * cfg.mat.thickness;
    Lr = max(cfg.geom.Lx_tot, cfg.geom.Ly_tot) / sqrt(12);     % rms moment arm
    lamC = max(eig(S.raw.C_hom));
    rows = [];  rhs = [];  chan = [];  nInfo = zeros(1,4);
    for ic = cases(:).'
        for iN = iNs(:).'
            T2 = S.T2(iN,ic);  T3 = S.T3(iN,ic);  T5 = S.T5(iN,ic);  T4 = S.T4(iN,ic);
            l2 = S.ell2(iN);  E = S.twoU_het(iN,ic);
            Sfe = squeeze(S.Sfe(iN,ic,:));  S2 = squeeze(S.S2(iN,ic,:));  S3 = squeeze(S.S3(iN,ic,:));
            Mfe = squeeze(S.Mfe(iN,ic,:));  M2 = squeeze(S.M2(iN,ic,:));  M3 = squeeze(S.M3(iN,ic,:));
            S4  = squeeze(S.S4(iN,ic,:));
            sigRef = sqrt(E*lamC/V) * V;

            % energy
            rows = [rows; [T2, T3, 2*T5, l2*T4] / E * sqrt(id.wE)];      %#ok<AGROW>
            rhs  = [rhs;  sqrt(id.wE)];                                   %#ok<AGROW>
            chan = [chan; 1];  nInfo(1) = nInfo(1) + 1;                   %#ok<AGROW>

            % stress
            nS = norm(Sfe);
            if nS >= id.tolInfo * sigRef
                w = sqrt(id.wS) / nS;
                rows = [rows; w*[zeros(3,1), S3, S2, zeros(3,1)]];        %#ok<AGROW>
                rhs  = [rhs;  w*Sfe];   chan = [chan; 2*ones(3,1)];       %#ok<AGROW>
                nInfo(2) = nInfo(2) + 1;
                if strcmpi(id.stressModel, 'both')
                    rows = [rows; w*[S2, zeros(3,1), S3, zeros(3,1)]];    %#ok<AGROW>
                    rhs  = [rhs;  w*Sfe];   chan = [chan; 3*ones(3,1)];   %#ok<AGROW>
                    nInfo(3) = nInfo(3) + 1;
                end
            end

            % double stress
            nM = norm(Mfe);
            if nM >= id.tolInfo * sigRef * Lr
                w = sqrt(id.wD) / nM;
                if strcmpi(id.doubleStressModel, 'mu')
                    blk = [zeros(6,1), zeros(6,1), zeros(6,1), l2*S4];
                else
                    blk = [zeros(6,1), M3, M2, l2*S4];
                end
                rows = [rows; w*blk];  rhs = [rhs; w*Mfe];  chan = [chan; 4*ones(6,1)]; %#ok<AGROW>
                nInfo(4) = nInfo(4) + 1;
            end
        end
    end
    wq  = [id.wE id.wS id.wS id.wD];
    use = wq(chan).' > 0;
    if id.fixAlphaM, free = [2 3 4];  am = id.alphaM;  else, free = [1 2 3 4];  am = NaN;  end
    Mf = rows(use, free);  bf = rhs(use);
    if id.fixAlphaM, bf = bf - rows(use,1)*am; end

    xu = Mf \ bf;
    a  = assemble(xu, free, am);
    method = 'unconstrained least squares';
    if ~feasible(a, id.pdMargin)
        [x, method] = constrainedSolve(Mf, bf, xu, free, am, id.pdMargin);
        a = assemble(x, free, am);
    end
    [ok, why] = feasible(a, id.pdMargin*0.999);

    res = rows*a(:) - rhs;
    rmsCh = nan(1,4);
    for q = 1:4
        m = use & chan == q;
        if any(m), rmsCh(q) = sqrt(mean((res(m)/sqrt(wq(q))).^2)); end
    end
    % channel-only estimates: each channel alone, on the columns it contains
    est = nan(4,4);
    for q = 1:4
        m = use & chan == q;
        if ~any(m), continue; end
        cols = find(any(rows(m,:) ~= 0, 1));
        if id.fixAlphaM, cols = setdiff(cols, 1); end
        bq = rhs(m);  if id.fixAlphaM, bq = bq - rows(m,1)*am; end
        if numel(cols) > 0 && nnz(m) >= numel(cols)
            est(q, cols) = (rows(m, cols) \ bq).';
        end
        if id.fixAlphaM, est(q,1) = am; end
    end
    r = struct('a', a, 'method', method, 'ok', ok, 'why', why, 'rmsCh', rmsCh, ...
               'sv', svd(Mf ./ max(abs(Mf), [], 1)), 'est', est, 'nInfo', nInfo);
end

function a = assemble(x, free, am)
    a = zeros(1,4);  a(free) = x(:).';
    if ~ismember(1, free), a(1) = am; end
    tiny = 1e-10*max(1, max(abs(a)));        % round-off negatives of bounded parameters
    if a(2) < 0 && a(2) > -tiny, a(2) = 0; end
    if a(4) < 0 && a(4) > -tiny, a(4) = 0; end
end

function [ok, why] = feasible(a, margin)
    why = {};
    if a(1) <= 0, why{end+1} = 'alpha_m <= 0'; end
    if a(2) <  0, why{end+1} = 'alpha_h < 0'; end
    if a(4) <  0, why{end+1} = 'alpha_l < 0'; end
    if a(1)*a(2) - a(3)^2 < margin*a(1)*a(2), why{end+1} = 'alpha_m*alpha_h - alpha_c^2 below margin'; end
    ok = isempty(why);  why = strjoin(why, ', ');
end

function [x, method] = constrainedSolve(M, b, x0, free, am, margin)
    H = M.'*M;  f = -M.'*b;
    nf = numel(free);
    lb = -inf(nf,1);  ub = inf(nf,1);
    ih = find(free == 2);  il = find(free == 4);  ic = find(free == 3);  im = find(free == 1);
    lb(ih) = 0;  lb(il) = 0;  if ~isempty(im), lb(im) = 1e-8; end
    if isempty(im), getm = @(x) am; else, getm = @(x) x(im); end
    xs = min(max(x0(:), lb + 1e-6), ub);             % feasible start
    xs(ih) = max(xs(ih), 1e-6);
    lim = sqrt(max((1 - margin)*getm(xs)*xs(ih), 0));
    xs(ic) = max(min(xs(ic), 0.99*lim), -0.99*lim);
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
            'SpecifyObjectiveGradient', true, 'OptimalityTolerance', 1e-12, ...
            'StepTolerance', 1e-14, 'MaxIterations', 2000);
        x = fmincon(@(x) qpObjective(x, H, f), xs, [], [], [], [], lb, ub, ...
            @(x) pdConstraint(x, ic, ih, im, am, margin), opts);
        method = 'fmincon (sqp), PD bound';
    elseif exist('sqp', 'file') == 2 || exist('sqp', 'builtin') == 5
        x = sqp(xs, {@(x) qpObjective(x, H, f), @(x) H*x(:) + f}, [], ...
                @(x) -pdConstraint(x, ic, ih, im, am, margin), lb, ub, 500, 1e-12);   % Octave
        method = 'sqp (Octave), PD bound';
    else
        x = x0;  method = 'unconstrained (no constrained solver available)';
    end
end

function [J, g] = qpObjective(x, H, f)
    % 0.5 x'Hx + f'x ; gradient only when requested (fmincon may ask for 1 output)
    x = x(:);
    J = 0.5*(x.'*H*x) + f.'*x;
    if nargout > 1, g = H*x + f; end
end

function [c, ceq] = pdConstraint(x, ic, ih, im, am, margin)
    % positive-definite local energy: am*ah - ac^2 >= margin*am*ah  <=>  c <= 0
    if isempty(im), m = am; else, m = x(im); end
    c = x(ic)^2 - (1 - margin)*m*x(ih);
    ceq = [];
end

% ============================================================
% ===================== trend alpha(N) =======================
% ============================================================

function T = fitTrends(N, aN, id)
    % alpha(N) = c0 + c1/N (+ c2/N^2 with >= 5 values of N), least squares
    N = N(:);  nK = numel(N);
    if nK >= 5, P = [0 1 2]; elseif nK >= 3, P = [0 1]; else, P = 0; end
    B = N .^ (-P);
    T.exponents = P;
    T.coef = nan(numel(P), 4);  T.rms = nan(1,4);  T.text = cell(1,4);
    for j = 1:4
        y = aN(:,j);
        if id.fixAlphaM && j == 1
            T.coef(:,j) = [id.alphaM; zeros(numel(P)-1,1)];  T.rms(j) = 0;
        else
            T.coef(:,j) = B \ y;
            T.rms(j) = sqrt(mean((B*T.coef(:,j) - y).^2));
        end
        parts = arrayfun(@(k) termText(T.coef(k,j), P(k)), 1:numel(P), 'UniformOutput', false);
        T.text{j} = strtrim(strjoin(parts, ' '));
    end
    fstr = {'c0', 'c1/N', 'c2/N^2'};
    T.formula = strjoin(fstr(1:numel(P)), ' + ');
    T.eval = @(Nq, j) (Nq(:) .^ (-P)) * T.coef(:,j);
end

function s = termText(c, p)
    if p == 0, s = sprintf('%.5g', c);
    elseif p == 1, s = sprintf('%+.5g/N', c);
    else, s = sprintf('%+.5g/N^%d', c, p);
    end
end

% ============================================================
% ======================== output ============================
% ============================================================

function writeAlphaCSV(fname, A)
    fid = fopen(fname, 'wt');
    if fid < 0, return; end
    fprintf(fid, 'N,alpha_m,alpha_h,alpha_c,alpha_l,rms_energy,rms_stress_psi,rms_stress_gamma,rms_double,PD_ok\n');
    if isfield(A, 'perN')
        for iN = 1:numel(A.Nlist)
            fprintf(fid, '%d', A.Nlist(iN));
            fprintf(fid, ',%.10g', A.perN.alphas(iN,:), A.perN.rmsPerChannel(iN,:));
            fprintf(fid, ',%d\n', A.perN.feasible(iN));
        end
    end
    if isfield(A, 'global')
        fprintf(fid, 'all');
        fprintf(fid, ',%.10g', A.global.alphas, A.global.rmsPerChannel);
        fprintf(fid, ',%d\n', A.global.feasible);
    end
    fclose(fid);
end

function plotAlphasVsN(A, id)
    N = A.Nlist(:);  aN = A.perN.alphas;
    Nf = linspace(min(N), max(N), 200).';
    lbl = {'\alpha_m', '\alpha_h', '\alpha_c', '\alpha_\ell'};
    figure('Color','w', 'Name', 'Alphas_vsN', 'Position', [40 60 1400 800]);
    onB = ~strcmp(A.perN.method, 'unconstrained least squares');
    for j = 1:4
        ax = subplot(2,2,j);  hold(ax,'on');  grid(ax,'on');
        h1 = plot(ax, N, aN(:,j), 'ko', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 7);
        h2 = plot(ax, Nf, A.trend.eval(Nf, j), 'b-', 'LineWidth', 1.4);
        hh = [h1 h2];  ll = {'identified at each N', ['trend ' A.trend.text{j}]};
        if any(onB)
            hc = plot(ax, N(onB), aN(onB,j), 'ms', 'MarkerSize', 11, 'LineWidth', 1.2);
            hh(end+1) = hc;  ll{end+1} = 'PD bound active'; %#ok<AGROW>
        end
        bad = ~A.perN.feasible;
        if any(bad)
            hb = plot(ax, N(bad), aN(bad,j), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
            hh(end+1) = hb;  ll{end+1} = 'PD violated'; %#ok<AGROW>
        end
        if isfield(A, 'global')
            hg = plot(ax, [min(N) max(N)], A.global.alphas(j)*[1 1], 'k--');
            hh(end+1) = hg;  ll{end+1} = 'one set for all N'; %#ok<AGROW>
        end
        set(ax, 'XTick', N);  xlim(ax, [min(N)-0.3, max(N)+0.3]);
        xlabel(ax, 'N (RUCs per side)');  ylabel(ax, lbl{j});
        ttl = sprintf('%s(N)', lbl{j});
        if id.fixAlphaM && j == 1, ttl = [ttl '  (fixed)']; end
        title(ax, ttl);
        legend(ax, hh, ll, 'Location', 'best');
    end
    sgtitle('Identified micromorphic parameters as functions of N');
end

function plotFit(S, aN, cases, figName, ttl, id, cfg)
    % model vs FE per channel; aN: nK x 4 alphas (one row per N); only the
    % groups used in the identification (informative) are shown
    Nk = S.Nlist(:);  nC = numel(cases);  cm = 0.85*hsv(nC);
    V  = cfg.geom.Lx_tot * cfg.geom.Ly_tot * cfg.mat.thickness;
    Lr = max(cfg.geom.Lx_tot, cfg.geom.Ly_tot) / sqrt(12);
    lamC = max(eig(S.raw.C_hom));
    figure('Color','w', 'Name', figName, 'Position', [30 30 1650 900]);
    axE = subplot(2,3,1);  hold(axE,'on');  grid(axE,'on');
    axS = subplot(2,3,2);  hold(axS,'on');  grid(axS,'on');
    axM = subplot(2,3,3);  hold(axM,'on');  grid(axM,'on');
    axSc = subplot(2,3,4); hold(axSc,'on'); grid(axSc,'on');
    axMc = subplot(2,3,5); hold(axMc,'on'); grid(axMc,'on');
    for t = 1:nC
        ic = cases(t);
        eE = nan(numel(Nk),1);  eS = eE;  eM = eE;
        for iN = 1:numel(Nk)
            a = aN(iN,:);
            Em = a(1)*S.T2(iN,ic) + a(2)*S.T3(iN,ic) + 2*a(3)*S.T5(iN,ic) + a(4)*S.ell2(iN)*S.T4(iN,ic);
            eE(iN) = Em/S.twoU_het(iN,ic) - 1;
            Sfe = squeeze(S.Sfe(iN,ic,:));  Mfe = squeeze(S.Mfe(iN,ic,:));
            Sm  = a(2)*squeeze(S.S3(iN,ic,:)) + a(3)*squeeze(S.S2(iN,ic,:));
            if strcmpi(id.doubleStressModel, 'mu')
                Mm = a(4)*S.ell2(iN)*squeeze(S.S4(iN,ic,:));
            else
                Mm = a(2)*squeeze(S.M3(iN,ic,:)) + a(3)*squeeze(S.M2(iN,ic,:)) + a(4)*S.ell2(iN)*squeeze(S.S4(iN,ic,:));
            end
            sigRef = sqrt(S.twoU_het(iN,ic)*lamC/V) * V;
            if norm(Sfe) >= id.tolInfo*sigRef
                eS(iN) = norm(Sm - Sfe)/norm(Sfe);
                plot(axSc, Sfe, Sm, '.', 'Color', cm(t,:), 'MarkerSize', 8);
            end
            if norm(Mfe) >= id.tolInfo*sigRef*Lr
                eM(iN) = norm(Mm - Mfe)/norm(Mfe);
                plot(axMc, Mfe, Mm, '.', 'Color', cm(t,:), 'MarkerSize', 8);
            end
        end
        plot(axE, Nk, eE, '-o', 'Color', cm(t,:), 'MarkerSize', 4);
        plot(axS, Nk, eS, '-o', 'Color', cm(t,:), 'MarkerSize', 4);
        plot(axM, Nk, eM, '-o', 'Color', cm(t,:), 'MarkerSize', 4);
    end
    title(axE, '(a) energy: W_{model}/2U_{het} - 1');
    title(axS, '(b) stress: |\int\sigma_{model} - \int\sigma_{FE}| / |\int\sigma_{FE}|');
    title(axM, '(c) double stress: |M_{model} - M_{FE}| / |M_{FE}|');
    for ax = [axE axS axM], set(ax, 'XTick', Nk); xlabel(ax, 'N'); end
    set(axS, 'YScale', 'log');  set(axM, 'YScale', 'log');
    for ax = [axSc axMc]
        lim = [min([get(ax,'XLim') get(ax,'YLim')]), max([get(ax,'XLim') get(ax,'YLim')])];
        plot(ax, lim, lim, 'k--');  xlabel(ax, 'FE');  ylabel(ax, 'model');  axis(ax, 'square');
    end
    title(axSc, '(d) \int\sigma components, all cases and N');
    title(axMc, '(e) \int\sigma_{ij} x_k components, all cases and N');
    sgtitle(ttl);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
