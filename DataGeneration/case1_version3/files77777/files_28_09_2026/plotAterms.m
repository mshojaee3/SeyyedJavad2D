function plotAterms(S, cfg)
% PLOTATERMS  One PNG per load case: A1..A19 versus N, raw (o) and smoothed (line).
%
%   plotAterms(S, cfg)          S from smoothData.m (or load('SmoothedData.mat'))
%
% Normalisation (per panel), to check the data after homogenization:
%   A1..A4   (psi psi)  :  A_j / AM_j        AM = same integral with psi -> eps_M
%                          e.g. A1 / int eps_M11^2 = C11 A1 / int eps_M11 C11 eps_M11
%   A5..A8   (gam gam)  :  A_j / AM_(j-4)    (size of gamma relative to eps_M)
%   A9..A13  (psi gam)  :  A9/AM1, A10/AM2, A11/AM2, A12/AM3, A13/AM4
%   A14..A19 (kap kap)  :  A_j / AM_j        AM = same integral with kappa -> grad eps_M
% psi -> eps_M means ratio 1 for A1..A4 (and A14..A19), 0 for A5..A13.
% The ratio is used only where the reference is physically relevant:
%   psi-type references: they carry >= 1% of the macroscopic energy,
%       c_ref |AM_ref| / 2U_hom >= 0.01  (c_ref = C11, 2|C12|, C22, 4 C33)
%   kappa references: ell^2 C11 AM_j / 2U_hom >= 1e-6 (not round-off) and
%       AM_j >= 1e-3 sum(AM_14..19)
% Otherwise (e.g. eps_M11 = 0 under H22, eps_M22 only an equilibrium
% correction under G1_22, grad eps_M = 0 under pure H loads) the panel shows the
% energy fraction  c_j A_j / 2U_hom  (c_j = weight of A_j in 2W with C_hom;
% kappa terms: ell^2 C11), and says so in its label.
%
% Files: <cfg.outDir>/A_terms/<case>_A.png   (cfg.outDir: region subfolder)

outDir = fullfile(cfg.outDir, 'A_terms');
if ~isfolder(outDir), mkdir(outDir); end
D    = S.raw;                        % region-selected raw data
keep = S.keep(:);
N    = S.Nlist(:);  Nr = D.Nlist(:);  Nr = Nr(keep);
nC   = size(S.A, 2);
C    = S.C_hom;
if ~isfield(D, 'AM')
    warning('plotAterms: no macroscopic references (AM) in the data -- rerun dataGeneration.');
    return;
end

ttl = {'A_{1} = \int\psi_{11}^2', 'A_{2} = \int\psi_{11}\psi_{22}', 'A_{3} = \int\psi_{22}^2', ...
       'A_{4} = \int\psi_{12}^2', 'A_{5} = \int\gamma_{11}^2', 'A_{6} = \int\gamma_{11}\gamma_{22}', ...
       'A_{7} = \int\gamma_{22}^2', 'A_{8} = \int\gamma_{12}^2', 'A_{9} = \int\psi_{11}\gamma_{11}', ...
       'A_{10} = \int\psi_{11}\gamma_{22}', 'A_{11} = \int\psi_{22}\gamma_{11}', ...
       'A_{12} = \int\psi_{22}\gamma_{22}', 'A_{13} = \int\psi_{12}\gamma_{12}', ...
       'A_{14} = \int\kappa_{111}^2', 'A_{15} = \int\kappa_{112}^2', 'A_{16} = \int\kappa_{122}^2', ...
       'A_{17} = \int\kappa_{222}^2', 'A_{18} = \int\kappa_{121}^2', 'A_{19} = \int\kappa_{221}^2'};
refIdx = [1 2 3 4, 1 2 3 4, 1 2 2 3 4, 14 15 16 17 18 19];     % AM used for A_j
refTxt = {'\int\epsilon_{M11}^2', '\int\epsilon_{M11}\epsilon_{M22}', '\int\epsilon_{M22}^2', ...
          '\int\epsilon_{M12}^2'};
cW = [C(1,1) 2*C(1,2) C(2,2) 4*C(3,3), C(1,1) 2*C(1,2) C(2,2) 4*C(3,3), ...
      2*C(1,1) 2*C(1,2) 2*C(1,2) 2*C(2,2) 8*C(3,3), C(1,1)*ones(1,6)];   % weights in 2W (C = C_hom)
ell2 = S.ell2(:);
reg  = 'full';  if isfield(S, 'region'), reg = S.region; end

for ic = 1:nC
    fh = figure('Color','w', 'Name', ['A_terms_' S.caseNames{ic}], 'Position', [20 20 1800 1050], ...
                'Visible', 'off');
    AMc = squeeze(D.AM(keep, ic, :));  if isvector(AMc), AMc = AMc(:).'; end
    U0 = S.twoU_hom(:, ic);
    cRef = [C(1,1) 2*abs(C(1,2)) C(2,2) 4*C(3,3)];          % weights of AM1..AM4 in 2W
    kapSum = sum(abs(AMc(:,14:19)), 2);
    for j = 1:19
        ax = subplot(4, 5, j);  hold(ax, 'on');  grid(ax, 'on');
        rawA = D.A(keep, ic, j);  smA = S.A(:, ic, j);
        den  = AMc(:, refIdx(j));
        if j <= 13
            useRef = all(cRef(refIdx(j)) * abs(den) ./ U0 >= 0.01);
        else
            useRef = all(ell2 * C(1,1) .* abs(den) ./ U0 >= 1e-6) && ...
                     all(abs(den) >= 1e-3 * kapSum);
        end
        if useRef
            yr = rawA ./ den;  ys = smA ./ den;
            if j <= 13
                lab = sprintf('A_{%d} / %s', j, refTxt{refIdx(j)});
            else
                lab = sprintf('A_{%d} / (same with \\nabla\\epsilon_M)', j);
            end
            ref = double(j <= 4 || j >= 14);         % psi -> eps_M: 1 (psi, kappa), 0 (gamma)
        else
            w = cW(j) * ones(size(N));
            if j >= 14, w = w .* ell2; end
            yr = w .* rawA ./ U0;  ys = w .* smA ./ U0;
            if j >= 14, lab = sprintf('\\ell^2 C_{11} A_{%d} / 2U_{hom}', j);
            else, lab = sprintf('c A_{%d} / 2U_{hom}', j); end
            ref = NaN;
        end
        plot(ax, Nr, yr, 'o', 'Color', [0.85 0.33 0.10], 'MarkerSize', 5);
        plot(ax, N,  ys, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4);
        if ~isnan(ref), plot(ax, [N(1) N(end)], [ref ref], 'k:'); end
        set(ax, 'XTick', N, 'FontSize', 7);  xlim(ax, [N(1)-0.2, N(end)+0.2]);
        title(ax, ttl{j}, 'FontSize', 8);  ylabel(ax, lab, 'FontSize', 7);
        if j > 14, xlabel(ax, 'N'); end
    end
    ax = subplot(4, 5, 20);  axis(ax, 'off');
    txt = {sprintf('%s   [%s]', strrep(S.caseNames{ic}, '_', '\_'), strrep(S.caseLabels{ic}, '_', '\_')), ...
           sprintf('region: %s', reg), '', 'o  raw data', '-   smoothed (used in the fit)', ...
           ':   value for \psi = \epsilon_M', '', ...
           'normalised by the same integral', 'of the macroscopic field \epsilon_M', ...
           '(energy fraction if that is 0)'};
    text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontSize', 9);
    sgtitle(fh, sprintf('A_1 ... A_{19} versus N  --  %s  [%s], region %s', ...
        strrep(S.caseNames{ic}, '_', '\_'), strrep(S.caseLabels{ic}, '_', '\_'), reg));
    print(fh, fullfile(outDir, sprintf('%s_A.png', S.caseNames{ic})), '-dpng', '-r130');
    close(fh);
end
fprintf('[plotAterms] %d PNG files written to %s\n', nC, outDir);
end
