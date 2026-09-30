function plotIterms(S, cfg)
% PLOTITERMS  One PNG per load case: I1..I24 versus N, raw (o) and smoothed (line).
%
%   plotIterms(S, cfg)          S from smoothData.m (or load('SmoothedData.mat'))
%
%   I1..I3   int psi_ij            I4..I6   int gamma_ij
%   I7..I12  int kappa_ijk         I13..I18 int psi_ij x_k      I19..I24 int gamma_ij x_k
%
% Normalisation (per panel): by the same integral of the macroscopic field eps_M,
%   I1..3 / IM1..3,   I4..6 / IM1..3,   I7..12 / IM7..12 (grad eps_M),
%   I13..18 / IM13..18,   I19..24 / IM13..18,
% so that psi -> eps_M gives 1 for the psi and kappa terms and 0 for the gamma
% terms. int gamma (I4..6) must be ~0 at every N (average-strain theorem of the
% KUBC: int eps = V sym(H) for every shift) -- an accuracy check of the data.
% The ratio is used only if the reference is relevant, |IM| >= 1e-2 * scale,
% with the scales  V eps_ref (I1..6),  V eps_ref / L_R (I7..12),  V L_r eps_ref
% (moments), eps_ref = sqrt(2U_hom / (C11 V)), L_R = sqrt(12) L_r; otherwise
% the panel shows I / scale and says so in its label.
%
% Files: <cfg.outDir>/I_terms/<case>_I.png   (cfg.outDir: region subfolder)

outDir = fullfile(cfg.outDir, 'I_terms');
if ~isfolder(outDir), mkdir(outDir); end
D    = S.raw;
keep = S.keep(:);
N    = S.Nlist(:);  Nr = D.Nlist(:);  Nr = Nr(keep);
nC   = size(S.I, 2);
C    = S.C_hom;
if ~isfield(D, 'IM')
    warning('plotIterms: no macroscopic references (IM) in the data -- rerun dataGeneration.');
    return;
end
V  = S.VolN(:);  Lr = S.LrN(:);  LR = sqrt(12)*Lr;

ttl = {'I_{1} = \int\psi_{11}', 'I_{2} = \int\psi_{22}', 'I_{3} = \int\psi_{12}', ...
       'I_{4} = \int\gamma_{11}', 'I_{5} = \int\gamma_{22}', 'I_{6} = \int\gamma_{12}', ...
       'I_{7} = \int\kappa_{111}', 'I_{8} = \int\kappa_{112}', 'I_{9} = \int\kappa_{221}', ...
       'I_{10} = \int\kappa_{222}', 'I_{11} = \int\kappa_{121}', 'I_{12} = \int\kappa_{122}', ...
       'I_{13} = \int\psi_{11}x_1', 'I_{14} = \int\psi_{11}x_2', 'I_{15} = \int\psi_{22}x_1', ...
       'I_{16} = \int\psi_{22}x_2', 'I_{17} = \int\psi_{12}x_1', 'I_{18} = \int\psi_{12}x_2', ...
       'I_{19} = \int\gamma_{11}x_1', 'I_{20} = \int\gamma_{11}x_2', 'I_{21} = \int\gamma_{22}x_1', ...
       'I_{22} = \int\gamma_{22}x_2', 'I_{23} = \int\gamma_{12}x_1', 'I_{24} = \int\gamma_{12}x_2'};
refIdx = [1 2 3, 1 2 3, 7:12, 13:18, 13:18];
refTxt = containers.Map('KeyType', 'double', 'ValueType', 'char');
names  = {'\epsilon_{M11}','\epsilon_{M22}','\epsilon_{M12}'};
for q = 1:3, refTxt(q) = ['\int' names{q}]; end
kn = {'\partial_1\epsilon_{M11}','\partial_2\epsilon_{M11}','\partial_1\epsilon_{M22}', ...
      '\partial_2\epsilon_{M22}','\partial_1\epsilon_{M12}','\partial_2\epsilon_{M12}'};
for q = 1:6, refTxt(6+q) = ['\int' kn{q}]; end
xn = {'x_1','x_2'};
for q = 1:6, refTxt(12+q) = sprintf('\\int%s%s', names{ceil(q/2)}, xn{2-mod(q,2)}); end
isGam = [false(1,3), true(1,3), false(1,12), true(1,6)];
reg = 'full';  if isfield(S, 'region'), reg = S.region; end

for ic = 1:nC
    fh = figure('Color','w', 'Name', ['I_terms_' S.caseNames{ic}], 'Position', [20 20 2000 1050], ...
                'Visible', 'off');
    IMc = squeeze(D.IM(keep, ic, :));  if isvector(IMc), IMc = IMc(:).'; end
    epsRef = sqrt(max(S.twoU_hom(:, ic), realmin) ./ (C(1,1) * V));
    for j = 1:24
        ax = subplot(4, 6, j);  hold(ax, 'on');  grid(ax, 'on');
        rawI = D.I(keep, ic, j);  smI = S.I(:, ic, j);
        if j <= 6,      sc = V .* epsRef;         scTxt = 'V\epsilon_{ref}';
        elseif j <= 12, sc = V .* epsRef ./ LR;   scTxt = 'V\epsilon_{ref}/L';
        else,           sc = V .* Lr .* epsRef;   scTxt = 'V L_r\epsilon_{ref}';
        end
        den = IMc(:, refIdx(j));
        useRef = all(abs(den) >= 1e-2 * sc);
        if useRef
            yr = rawI ./ den;  ys = smI ./ den;
            lab = sprintf('I_{%d} / %s', j, refTxt(refIdx(j)));
            ref = double(~isGam(j));
        else
            yr = rawI ./ sc;  ys = smI ./ sc;
            lab = sprintf('I_{%d} / %s', j, scTxt);
            ref = NaN;  if isGam(j), ref = 0; end
        end
        plot(ax, Nr, yr, 'o', 'Color', [0.85 0.33 0.10], 'MarkerSize', 5);
        plot(ax, N,  ys, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4);
        if ~isnan(ref), plot(ax, [N(1) N(end)], [ref ref], 'k:'); end
        set(ax, 'XTick', N, 'FontSize', 7);  xlim(ax, [N(1)-0.2, N(end)+0.2]);
        title(ax, ttl{j}, 'FontSize', 8);  ylabel(ax, lab, 'FontSize', 7);
        if j > 18, xlabel(ax, 'N'); end
    end
    sgtitle(fh, sprintf(['I_1 ... I_{24} versus N  --  %s  [%s], region %s     ' ...
        '(o raw, - smoothed, : value for \\psi = \\epsilon_M)'], ...
        strrep(S.caseNames{ic}, '_', '\_'), strrep(S.caseLabels{ic}, '_', '\_'), reg));
    print(fh, fullfile(outDir, sprintf('%s_I.png', S.caseNames{ic})), '-dpng', '-r130');
    close(fh);
end
fprintf('[plotIterms] %d PNG files written to %s\n', nC, outDir);
end
