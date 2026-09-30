% ============================================================
%  smoothEnergyTerms_vsN.m
%
%  Post-processing of EnergyTerms_vsN.mat (written by
%  MicroMacroStrain_2D_LE22.m): removes the numerical fluctuation of the
%  energy terms over N so that every term is monotone in N, for use in
%  the alpha identification. No FE solve is repeated.
%
%  What is noise and what is not (per load case, C = C_hom):
%   T1 = int eps_M' C eps_M   Stage-1 field, independent of N  -> constant
%   T3 = int gamma' C gamma   positive, decays with N; small scatter
%                             (few %) from the finite shift set
%                             -> monotone fit, see below
%   T5 = int psi' C gamma     continuum: T3 + T5 = int eps_M' C gamma = 0
%                             (eps_M equilibrated, gamma from a displacement
%                             that vanishes on the KUBC boundary). The raw T5
%                             carries the discrete residual of that identity,
%                             which is of the order of T3 itself
%                             -> T5 := -T3  (exact, no fit)
%   T2 = int psi' C psi       T2 - T1 is a difference of two large grid
%                             integrals, dominated by quadrature error
%                             -> T2 := T1 + T3 (exact, no fit)
%   T4 = int kappa' A kappa   tends to the macro curvature (N-independent);
%                             +/-1% scatter -> non-increasing fit c0+c1/N+c2/N^2
%   2U_het                    heterogeneous energy, excess over 2U_hom decays
%                             -> 2U_hom + monotone fit of (2U_het - 2U_hom)
%
%  Monotone fit of a decaying quantity y(N):
%       y(N) = sum_k c_k N^(-p_k),  c_k >= 0,  p = decayExponents
%  (non-negative combination of decreasing functions => strictly
%  decreasing, convex, -> 0), relative least squares via lsqnonneg.
%
%  The smoothed terms satisfy T1 = T2 + 2 T5 + T3 and T3 + T5 = 0 exactly.
%
%  N = 1 is a different regime: the single cell is entirely inside the
%  KUBC boundary layer, and for some load cases (e.g. G1_11) T3 rises from
%  N = 1 to N = 2 -- a systematic effect, not noise. With NminFit = 2
%  (default) N = 1 is excluded from the fit and from the smoothed output.
% ============================================================

clear; clc;

%% ============ USER INPUT ============

inFile   = fullfile(pwd, 'MicroMacroStrain_out', 'EnergyTerms_vsN.mat');
outDir   = fileparts(inFile);

NminFit        = 2;              % fit and keep only N >= NminFit
decayExponents = [0.5 1 1.5 2];  % basis N^-p for T3 and 2U_het - 2U_hom
T4Exponents    = [0 1 2];        % basis N^-p for T4 (constant + decay)

plotCases      = [];             % [] -> all cases
plotPerCase    = true;           % 3-panel raw-vs-smoothed figure per case
SAVE_FIGS      = false;

%% ============ LOAD ============

VN = load(inFile);
Nl = VN.Nlist(:);
nN = numel(Nl);
[~, nC] = size(VN.T1);
if ~isfield(VN, 'caseLabels'), VN.caseLabels = VN.caseNames; end
if ~isfield(VN, 'twoU_hom'),   VN.twoU_hom   = VN.T1;        end
if isempty(plotCases), plotCases = 1:nC; end

keep = Nl >= NminFit;
if nnz(keep) < 3
    error('Need at least 3 N values >= NminFit (have %d).', nnz(keep));
end
ell2 = VN.ell2T4(:,1) ./ VN.T4(:,1);          % (L/N)^2 for each N

%% ============ SMOOTH ============

S = struct();
flds = {'T1','T2','T3','T4','T5','ell2T4','twoU_het','twoU_hom'};
for f = flds, S.(f{1}) = nan(nN, nC); end
S.coefT3 = nan(numel(decayExponents), nC);
S.coefDU = nan(numel(decayExponents), nC);
S.coefT4 = nan(numel(T4Exponents), nC);
Q = struct('rmsT3', nan(1,nC), 'rmsDU', nan(1,nC), 'rmsT4', nan(1,nC), ...
           'T1spread', nan(1,nC), 'orthRaw', nan(1,nC), 'T2mT1Raw', nan(1,nC));

for ic = 1:nC
    % T1: constant (Stage-1 field does not depend on N)
    T1c = mean(VN.T1(keep,ic));
    Q.T1spread(ic) = max(abs(VN.T1(keep,ic) - T1c)) / T1c;

    % T3: monotone decaying fit
    [T3s, S.coefT3(:,ic), Q.rmsT3(ic)] = monoFit(Nl, VN.T3(:,ic), keep, decayExponents);

    % T4: constant + decay (non-increasing)
    [T4s, S.coefT4(:,ic), Q.rmsT4(ic)] = monoFit(Nl, VN.T4(:,ic), keep, T4Exponents);

    % heterogeneous energy: exact 2U_hom + monotone excess
    U0 = mean(VN.twoU_hom(keep,ic));
    [DUs, S.coefDU(:,ic), Q.rmsDU(ic)] = monoFit(Nl, VN.twoU_het(:,ic) - U0, keep, decayExponents);

    S.T1(:,ic)       = T1c;
    S.T3(:,ic)       = T3s;
    S.T5(:,ic)       = -T3s;          % T3 + T5 = 0
    S.T2(:,ic)       = T1c + T3s;     % T1 = T2 + 2 T5 + T3
    S.T4(:,ic)       = T4s;
    S.ell2T4(:,ic)   = ell2 .* T4s;
    S.twoU_hom(:,ic) = U0;
    S.twoU_het(:,ic) = U0 + DUs;

    % size of the raw inconsistencies that the identities remove
    Q.orthRaw(ic)  = max(abs(VN.T3(keep,ic) + VN.T5(keep,ic)) ./ VN.T3(keep,ic));
    Q.T2mT1Raw(ic) = max(abs(VN.T2(keep,ic) - VN.T1(keep,ic) - VN.T3(keep,ic)) ./ VN.T3(keep,ic));
end
for f = flds, S.(f{1})(~keep,:) = NaN; end   % N < NminFit not used

% ---- checks: monotonicity and identities ----
dT3 = diff(S.T3(keep,:));  dT4 = diff(S.T4(keep,:));  dU = diff(S.twoU_het(keep,:));
monoOK = all(dT3(:) < 0) && all(dT4(:) <= 0) && all(dU(:) < 0);
idErr  = max(max(abs(S.T1(keep,:) - S.T2(keep,:) - 2*S.T5(keep,:) - S.T3(keep,:)) ./ S.T1(keep,:)));

%% ============ REPORT ============

fprintf('Smoothing of %d load cases, N = %s used (NminFit = %d)\n', nC, mat2str(Nl(keep).'), NminFit);
fprintf('  T3, 2U_het-2U_hom: sum c_k N^-p, p = %s, c_k >= 0\n', mat2str(decayExponents));
fprintf('  T4               : sum c_k N^-p, p = %s, c_k >= 0\n', mat2str(T4Exponents));
fprintf('  T5 = -T3,  T2 = T1 + T3  (continuum identities)\n\n');
fprintf('  %-9s %-28s %8s %8s %8s %9s | %10s %11s\n', 'Case', 'label', 'rms T3', ...
    'rms dU', 'rms T4', 'T1 spread', 'raw|T3+T5|', 'raw|T2-T1-T3|');
fprintf('  %-9s %-28s %8s %8s %8s %9s | %10s %11s\n', '', '', '(rel)', '(rel)', '(rel)', '(rel)', '/T3 max', '/T3 max');
for ic = 1:nC
    fprintf('  %-9s %-28s %8.3f %8.3f %8.4f %9.1e | %10.2f %11.2f\n', VN.caseNames{ic}, ...
        VN.caseLabels{ic}, Q.rmsT3(ic), Q.rmsDU(ic), Q.rmsT4(ic), Q.T1spread(ic), ...
        Q.orthRaw(ic), Q.T2mT1Raw(ic));
end
fprintf('\n  monotone (T3, 2U_het strictly decreasing; T4 non-increasing): %s\n', mat2str(monoOK));
fprintf('  identity T1 = T2 + 2T5 + T3, max relative error: %.1e\n', idErr);

% ---- save ----
S.Nlist = Nl.';  S.keep = keep.';  S.NminFit = NminFit;
S.caseNames = VN.caseNames;  S.caseLabels = VN.caseLabels;
if isfield(VN, 'rawVecs'), S.rawVecs = VN.rawVecs; end
S.decayExponents = decayExponents;  S.T4Exponents = T4Exponents;
S.fitQuality = Q;  S.raw = VN;
save(fullfile(outDir, 'EnergyTerms_vsN_smoothed.mat'), '-struct', 'S');

fid = fopen(fullfile(outDir, 'EnergyTerms_vsN_smoothed.csv'), 'wt');
if fid > 0
    fprintf(fid, 'Case,LoadCase,Label,N,T1,T2,T3,T4,T5,ell2_T4,twoU_hom,twoU_het\n');
    for ic = 1:nC
        for iN = find(keep).'
            fprintf(fid, '%d,%s,"%s",%d', ic, VN.caseNames{ic}, VN.caseLabels{ic}, Nl(iN));
            fprintf(fid, ',%.10e', S.T1(iN,ic), S.T2(iN,ic), S.T3(iN,ic), S.T4(iN,ic), ...
                S.T5(iN,ic), S.ell2T4(iN,ic), S.twoU_hom(iN,ic), S.twoU_het(iN,ic));
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end
fprintf('  saved EnergyTerms_vsN_smoothed.mat / .csv in %s\n', outDir);

%% ============ PLOTS ============

fg = plotGridSmoothed(Nl, VN, S, plotCases);
if SAVE_FIGS, print(fg, fullfile(outDir, 'AllCases_ratios_vsN_smoothed.png'), '-dpng', '-r150'); end
if plotPerCase
    for ic = plotCases(:).'
        fc = plotCaseSmoothed(Nl, VN, S, ic);
        if SAVE_FIGS
            print(fc, fullfile(outDir, sprintf('%s_EnergyTerms_vsN_smoothed.png', VN.caseNames{ic})), ...
                '-dpng', '-r150');
        end
    end
end


% ============================================================
% ======================= LOCAL FUNCTIONS =====================
% ============================================================

function [ys, c, rmsRel] = monoFit(N, y, keep, P)
    % y(N) ~ sum_k c_k N^-P(k), c_k >= 0; relative least squares on keep
    N = N(:);  y = y(:);
    A = N .^ (-P(:).');
    w = 1 ./ abs(y(keep));
    if any(~isfinite(w)) || any(y(keep) <= 0)
        w = ones(nnz(keep),1) / max(abs(y(keep)));   % fall back to absolute scaling
    end
    c  = lsqnonneg(A(keep,:) .* w, y(keep) .* w);
    ys = A * c;
    rmsRel = sqrt(mean(((ys(keep) - y(keep)) ./ y(keep)).^2));
end

function fh = plotGridSmoothed(Nl, VN, S, cases)
    % one log-log tile per case: raw (markers) vs smoothed (lines), / T1
    K = numel(cases);  nc = ceil(sqrt(K));  nr = ceil(K/nc);
    col = [0.929 0.694 0.125; 0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933];
    k = S.keep(:);
    fh = figure('Color','w', 'Name', 'AllCases_ratios_vsN_smoothed', 'Position', [20 20 1800 1150]);
    allv = [];
    for ic = cases(:).'
        T1 = S.T1(k,ic);
        allv = [allv; VN.T3(k,ic)./T1; VN.ell2T4(k,ic)./T1; ...
                (VN.twoU_het(k,ic) - S.twoU_hom(k,ic))./T1]; %#ok<AGROW>
    end
    allv = allv(allv > 0);
    yl = [10^floor(log10(min(allv))), 10^ceil(log10(max(allv)))];
    for t = 1:K
        ic = cases(t);
        ax = subplot(nr, nc, t);  hold(ax,'on');  grid(ax,'on');
        T1r = VN.T1(:,ic);  T1s = S.T1(k,ic);
        rawv = {VN.T3(:,ic)./T1r, -VN.T5(:,ic)./T1r, VN.ell2T4(:,ic)./T1r, ...
                (VN.twoU_het(:,ic) - VN.twoU_hom(:,ic))./T1r};
        smv  = {S.T3(k,ic)./T1s, -S.T5(k,ic)./T1s, S.ell2T4(k,ic)./T1s, ...
                (S.twoU_het(k,ic) - S.twoU_hom(k,ic))./T1s};
        mk = {'^','v','d','p'};
        h = [];
        for q = [1 3 4]                       % smoothed -T5 coincides with T3
            h(end+1) = plot(ax, Nl(k), smv{q}, '-', 'Color', col(q,:), 'LineWidth', 1.6); %#ok<AGROW>
        end
        for q = 1:4
            v = rawv{q};  v(v <= 0) = NaN;
            plot(ax, Nl(k), v(k), mk{q}, 'Color', col(q,:), 'MarkerSize', 4);
            plot(ax, Nl(~k), v(~k), mk{q}, 'Color', [0.6 0.6 0.6], 'MarkerSize', 4);
        end
        set(ax, 'XScale','log', 'YScale','log', 'XTick', Nl, 'FontSize', 7);
        xlim(ax, [Nl(1)*0.9, Nl(end)*1.1]);  ylim(ax, yl);
        title(ax, sprintf('%d: %s', ic, strrep(VN.caseLabels{ic}, '_', '\_')), 'FontSize', 8);
        if t > (nr-1)*nc, xlabel(ax, 'N'); end
        if mod(t-1, nc) == 0, ylabel(ax, 'ratio to T1'); end
        if t == 1
            legend(ax, h, {'T3/T1 = -T5/T1', '(L/N)^2T4/T1', '(2U_{het}-2U_{hom})/T1'}, ...
                'Location', 'southwest', 'FontSize', 6);
        end
    end
    sgtitle(fh, sprintf(['Smoothed (lines) vs raw (markers: \\Delta T3, \\nabla -T5, ' ...
        '\\diamond (L/N)^2T4, \\star 2U excess; grey = N < %d, not used)'], S.NminFit));
end

function fh = plotCaseSmoothed(Nl, VN, S, ic)
    k = S.keep(:);
    col = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125; ...
           0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933];
    fh = figure('Color','w', 'Name', ['Smoothed_' VN.caseNames{ic}], 'Position', [40 80 1700 520]);
    T1s = S.T1(k,ic);

    % (a) terms
    ax = subplot(1,3,1);  hold(ax,'on');  grid(ax,'on');
    nm = {'T1','T2','T3','T5','ell2T4','twoU_het'};
    lb = {'T1','T2','T3','T5','(L/N)^2T4','2U_{het}'};
    h = [];
    for q = 1:6
        h(end+1) = plot(ax, Nl(k), S.(nm{q})(k,ic), '-', 'Color', col(q,:), 'LineWidth', 1.6); %#ok<AGROW>
        plot(ax, Nl, VN.(nm{q})(:,ic), 'o', 'Color', col(q,:), 'MarkerSize', 4);
    end
    legend(ax, h, lb, 'Location', 'best');
    xlabel(ax, 'N');  ylabel(ax, 'domain integral');  set(ax, 'XTick', Nl);
    title(ax, '(a) energy terms: raw (o), smoothed (line)');

    % (b) -> 1
    ax = subplot(1,3,2);  hold(ax,'on');  grid(ax,'on');
    h1 = plot(ax, Nl(k), S.T2(k,ic)./T1s, '-', 'Color', col(2,:), 'LineWidth', 1.6);
    plot(ax, Nl, VN.T2(:,ic)./VN.T1(:,ic), 'o', 'Color', col(2,:), 'MarkerSize', 4);
    h2 = plot(ax, Nl(k), S.twoU_het(k,ic)./T1s, '-', 'Color', col(6,:), 'LineWidth', 1.6);
    plot(ax, Nl, VN.twoU_het(:,ic)./VN.T1(:,ic), 'o', 'Color', col(6,:), 'MarkerSize', 4);
    plot(ax, [Nl(1) Nl(end)], [1 1], 'k:');
    legend(ax, [h1 h2], {'T2/T1', '2U_{het}/T1'}, 'Location', 'best');
    xlabel(ax, 'N');  ylabel(ax, 'ratio to T1');  set(ax, 'XTick', Nl);
    title(ax, '(b) should tend to 1');

    % (c) -> 0, log-log
    ax = subplot(1,3,3);  hold(ax,'on');  grid(ax,'on');
    sm = {S.T3(k,ic)./T1s, S.ell2T4(k,ic)./T1s, (S.twoU_het(k,ic) - S.twoU_hom(k,ic))./T1s};
    rw = {VN.T3(:,ic)./VN.T1(:,ic), VN.ell2T4(:,ic)./VN.T1(:,ic), ...
          (VN.twoU_het(:,ic) - VN.twoU_hom(:,ic))./VN.T1(:,ic)};
    cc = col([3 5 6],:);
    h = [];
    for q = 1:3
        h(end+1) = plot(ax, Nl(k), sm{q}, '-', 'Color', cc(q,:), 'LineWidth', 1.6); %#ok<AGROW>
        v = rw{q};  v(v <= 0) = NaN;
        plot(ax, Nl, v, 'o', 'Color', cc(q,:), 'MarkerSize', 4);
    end
    v = -VN.T5(:,ic)./VN.T1(:,ic);  v(v <= 0) = NaN;
    h(end+1) = plot(ax, Nl, v, 'v', 'Color', col(4,:), 'MarkerSize', 4);
    set(ax, 'XScale','log', 'YScale','log', 'XTick', Nl);
    xlim(ax, [Nl(1)*0.9, Nl(end)*1.1]);
    legend(ax, h, {'T3/T1 (= -T5/T1 smoothed)', '(L/N)^2T4/T1', '(2U_{het}-2U_{hom})/T1', ...
        'raw -T5/T1'}, 'Location', 'best');
    xlabel(ax, 'N');  ylabel(ax, 'ratio to T1');
    title(ax, '(c) should vanish (log-log)');

    sgtitle(fh, sprintf('Raw vs smoothed  --  %s  [%s]', strrep(VN.caseNames{ic}, '_', '\_'), ...
        strrep(VN.caseLabels{ic}, '_', '\_')));
end