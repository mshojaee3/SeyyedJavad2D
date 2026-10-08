function save_figures(P1, H, Dident, par1, res)
% =========================================================================
%  PNG FIGURES of the whole run (called by main_micromorphic after all stages; nothing here changes the saved data)
%      homogenization/figures/   displacement, strain and stress of the three periodic unit-strain problems
%      case1/figures/            identified parameters theta(N);  residuals;  fullscale/N<N>/: displacement and strain/stress (bin means before and
%                                after smoothing) of the heterogeneous plates and of the homogeneous plate
%      case2/figures/            fullscale/ (same as case 1),  micromorphic/ (displacement, eps_macro, psi, gamma, sigma),
%                                comparison/ (full-scale vs micromorphic vs homogeneous: displacement, psi, sigma, errors),  error_summary.png
%  A failing figure only gives a warning; the data are already saved.
% =========================================================================
dpi = P1.out.pngDPI;
homDir = P1.out.homFolder;
dir1 = P1.out.folder;

safe(@() fig_homogenization(H, fullfile(homDir, 'figures'), dpi), 'homogenization');
safe(@() fig_identification(P1, par1, fullfile(dir1, 'figures'), dpi), 'identification');
safe(@() fig_fullscale(P1, Dident, 'ident', fullfile(dir1, 'figures', 'fullscale'), dpi), 'full-scale, case 1');

for ic = 1:numel(res)
    r = res{ic};
    d2 = r.P2.out.folder;
    safe(@() fig_fullscale(r.P2, r.Dtest, 'test', fullfile(d2, 'figures', 'fullscale'), dpi), 'full-scale, case 2');
    safe(@() fig_micromorphic(r.P2, r.M, fullfile(d2, 'figures', 'micromorphic'), dpi), 'micromorphic');
    safe(@() fig_comparison(r.P2, r.Dtest, r.M, fullfile(d2, 'figures', 'comparison'), dpi), 'comparison fields');
    safe(@() fig_error_summary(r.C, fullfile(d2, 'figures'), dpi), 'error summary');
end
fprintf('  figures saved.\n');
end


function safe(fcn, label)
try
    fcn();
catch err
    warning('save_figures:failed', 'Figures "%s" failed: %s', label, err.message);
end
end


%% ========================================================================
%  FIGURE FUNCTIONS
% ========================================================================

function fig_homogenization(H, outDir, dpi)
info = H.info;  nodes = info.nodes;  elems = info.elems6;
eN = {'\epsilon_{11}', '\epsilon_{22}', '\gamma_{12}'};
sN = {'\sigma_{11}', '\sigma_{22}', '\sigma_{12}'};
for k = 1:3
    u = info.u{k};  f = info.fld{k};
    fig = new_fig([50 50 1600 760]);
    subplot(2, 4, 1);  mesh_panel(gca, nodes, elems, u(1:2:end), 'node', 'u_x', []);
    subplot(2, 4, 5);  mesh_panel(gca, nodes, elems, u(2:2:end), 'node', 'u_y', []);
    for c = 1:3
        subplot(2, 4, 1 + c);  mesh_panel(gca, nodes, elems, mean(f.eps(:, :, c), 2), 'elem', eN{c}, []);
        subplot(2, 4, 5 + c);  mesh_panel(gca, nodes, elems, mean(f.sig(:, :, c), 2), 'elem', sN{c}, []);
    end
    add_title(fig, sprintf('Homogenization of the periodic unit cell, unit macro strain %s', eN{k}));
    save_png(fig, fullfile(outDir, sprintf('homogenization_unit_strain_%d.png', k)), dpi);
end
end


function fig_identification(P, par, outDir, dpi)
Chom = P.material.Chom;
hom3 = [(Chom(1,1) + Chom(2,2))/2; Chom(1,2); Chom(3,3)];
NL = P.ident.NList;
Nmax = max(4*NL(end), 10);  Nv = linspace(1, Nmax, 200);
mats = {par.c_psi, par.c_gamma, par.c_couple};
mname = {'C_\psi', 'C_\gamma', 'C_{couple}'};
cname = {'a (C_{11} = C_{22})', 'b (C_{12})', 'd (C_{33})'};
fig = new_fig([50 50 1400 900]);
for m = 1:3
    for q = 1:3
        ax = subplot(3, 3, (m-1)*3 + q);  hold(ax, 'on');
        c = mats{m}(q, :);
        plot(ax, Nv, c(1) + c(2)./Nv + c(3)./Nv.^2, 'LineWidth', 1.8);
        plot(ax, [1 Nmax], hom3(q)*[1 1], '--k');
        yl = get(ax, 'YLim');
        plot(ax, NL(end)*[1 1], yl, ':', 'Color', [0.5 0.5 0.5]);
        set(ax, 'YLim', yl);
        grid(ax, 'on');  box(ax, 'on');
        title(ax, [mname{m} ':  ' cname{q}], 'FontSize', 9);
        if m == 3, xlabel(ax, 'N'); end
        if m == 1 && q == 1, legend(ax, {'identified \theta(N)', 'C_{hom}', 'last identified N'}, 'Location', 'best'); end
    end
end
add_title(fig, 'Identified micromorphic stiffnesses  \theta(N) = c_0 + c_1/N + c_2/N^2');
save_png(fig, fullfile(outDir, 'identified_parameters.png'), dpi);

if isfield(par, 'relRes') && isfield(par, 'groups')
    fig = new_fig([50 50 900 520]);
    ax = axes('Parent', fig);
    bar(ax, NL(:), par.relRes);
    grid(ax, 'on');  xlabel(ax, 'N');  ylabel(ax, 'relative residual');
    legend(ax, par.groups, 'Location', 'best');
    title(ax, 'Residuals of the theory per group (identified parameters)', 'FontSize', 10);
    save_png(fig, fullfile(outDir, 'identification_residuals.png'), dpi);
end
end


function fig_fullscale(P, D, which, outDir, dpi)
pl = P.(which);  Lx = pl.Lx;  Ly = pl.Ly;  loads = pl.loads;  Ns = pl.NList;
hm = D.hom.mesh;
for l = 1:numel(loads)
    % homogeneous plate (independent of N)
    ux = hm.U(1:2:end, l);  uy = hm.U(2:2:end, l);
    fig = new_fig([50 50 1300 420]);
    subplot(1, 2, 1);  mesh_panel(gca, hm.nodes, hm.elems, ux, 'node', 'u_x  homogeneous plate', []);
    subplot(1, 2, 2);  mesh_panel(gca, hm.nodes, hm.elems, uy, 'node', 'u_y  homogeneous plate', []);
    add_title(fig, sprintf('Homogeneous plate, load: %s', loads{l}));
    save_png(fig, fullfile(outDir, 'homogeneous', sprintf('%s_displacement.png', loads{l})), dpi);
    fields_figure(D.hom.epsBin(:,:,:,l), D.hom.epsSmooth(:,:,:,l), D.hom.sigBin(:,:,:,l), D.hom.sigSmooth(:,:,:,l), Lx, Ly, ...
        sprintf('Homogeneous plate, load: %s', loads{l}), fullfile(outDir, 'homogeneous', sprintf('%s_fields.png', loads{l})), dpi);
end
for iN = 1:numel(Ns)
    h = D.het(iN);  N = Ns(iN);  hmesh = h.mesh;
    for l = 1:numel(loads)
        ux1 = hmesh.U(1:2:end, l);  uy1 = hmesh.U(2:2:end, l);
        uxh = hm.U(1:2:end, l);     uyh = hm.U(2:2:end, l);
        limx = vlim([ux1; uxh]);  limy = vlim([uy1; uyh]);
        fig = new_fig([50 50 1500 620]);
        subplot(2, 2, 1);  mesh_panel(gca, hmesh.nodes, hmesh.elems, ux1, 'node', 'u_x  heterogeneous (full-scale)', limx);
        subplot(2, 2, 2);  mesh_panel(gca, hmesh.nodes, hmesh.elems, uy1, 'node', 'u_y  heterogeneous (full-scale)', limy);
        subplot(2, 2, 3);  mesh_panel(gca, hm.nodes, hm.elems, uxh, 'node', 'u_x  homogeneous', limx);
        subplot(2, 2, 4);  mesh_panel(gca, hm.nodes, hm.elems, uyh, 'node', 'u_y  homogeneous', limy);
        add_title(fig, sprintf('%s: N = %d, load: %s (centred inclusions)', which, N, loads{l}));
        dirN = fullfile(outDir, sprintf('N%d', N));
        save_png(fig, fullfile(dirN, sprintf('%s_displacement.png', loads{l})), dpi);
        fields_figure(h.epsBin(:,:,:,l), h.epsSmooth(:,:,:,l), h.sigBin(:,:,:,l), h.sigSmooth(:,:,:,l), Lx, Ly, ...
            sprintf('%s: N = %d, load: %s, bins of 1/%d cell (psiMode: %s)', which, N, loads{l}, P.fe.NB, pl.psiMode), ...
            fullfile(dirN, sprintf('%s_fields.png', loads{l})), dpi);
    end
end
end


function fields_figure(Er, Es, Sr, Ss, Lx, Ly, ttl, file, dpi)
% strain and stress: bin means (raw, before smoothing) and after the Legendre smoothing.   Inputs: nby x nbx x 3
eN = {'\epsilon_{11}', '\epsilon_{22}', '\gamma_{12}'};
sN = {'\sigma_{11}', '\sigma_{22}', '\sigma_{12}'};
fig = new_fig([50 50 1700 760]);
for c = 1:3
    a = Er(:, :, c);  b = Es(:, :, c);  le = vlim([a(:); b(:)]);
    subplot(3, 4, (c-1)*4 + 1);  grid_panel(gca, a, Lx, Ly, [eN{c} '  bin means (raw)'], le);
    subplot(3, 4, (c-1)*4 + 2);  grid_panel(gca, b, Lx, Ly, [eN{c} '  smoothed'], le);
    a = Sr(:, :, c);  b = Ss(:, :, c);  ls = vlim([a(:); b(:)]);
    subplot(3, 4, (c-1)*4 + 3);  grid_panel(gca, a, Lx, Ly, [sN{c} '  bin means (raw)'], ls);
    subplot(3, 4, (c-1)*4 + 4);  grid_panel(gca, b, Lx, Ly, [sN{c} '  smoothed'], ls);
end
add_title(fig, ttl);
save_png(fig, file, dpi);
end


function fig_micromorphic(P, M, outDir, dpi)
T = P.test;  loads = T.loads;  Ns = T.NList;
nodes = M.mesh.nodes;  elems = M.mesh.elems;
eN = {'\epsilon_{11}', '\epsilon_{22}', '\gamma_{12}'};
sN = {'\sigma_{11}', '\sigma_{22}', '\sigma_{12}'};
pN = {'\psi_{11}', '\psi_{22}', '\psi_{12}'};
gN = {'\gamma_{11}', '\gamma_{22}', '\gamma_{12}'};
for iN = 1:numel(Ns)
    raw = M.N(iN).raw;  N = Ns(iN);
    for l = 1:numel(loads)
        U = raw.U(:, l);
        fig = new_fig([50 50 1500 420]);
        subplot(1, 2, 1);  mesh_panel(gca, nodes, elems, U(1:5:end), 'node', 'u_x  micromorphic', []);
        subplot(1, 2, 2);  mesh_panel(gca, nodes, elems, U(2:5:end), 'node', 'u_y  micromorphic', []);
        add_title(fig, sprintf('Micromorphic solution: N = %d, load: %s', N, loads{l}));
        dirN = fullfile(outDir, sprintf('N%d', N));
        save_png(fig, fullfile(dirN, sprintf('%s_displacement.png', loads{l})), dpi);

        fig = new_fig([50 50 1700 760]);
        for c = 1:3
            ee = elem_mean(raw.Wg, raw.Eps(:, :, c, l));
            pp = elem_mean(raw.Wg, raw.Psi(:, :, c, l));
            ss = elem_mean(raw.Wg, raw.Sig(:, :, c, l));
            gg = ee - pp;
            subplot(3, 4, (c-1)*4 + 1);  mesh_panel(gca, nodes, elems, ee, 'elem', [eN{c} '  macro strain sym grad u'], []);
            subplot(3, 4, (c-1)*4 + 2);  mesh_panel(gca, nodes, elems, pp, 'elem', pN{c}, []);
            subplot(3, 4, (c-1)*4 + 3);  mesh_panel(gca, nodes, elems, gg, 'elem', [gN{c} ' = \epsilon - \psi'], []);
            subplot(3, 4, (c-1)*4 + 4);  mesh_panel(gca, nodes, elems, ss, 'elem', sN{c}, []);
        end
        add_title(fig, sprintf('Micromorphic fields (element means): N = %d, load: %s', N, loads{l}));
        save_png(fig, fullfile(dirN, sprintf('%s_fields.png', loads{l})), dpi);
    end
end
end


function fig_comparison(P, D, M, outDir, dpi)
T = P.test;  Lx = T.Lx;  Ly = T.Ly;  loads = T.loads;  Ns = T.NList;
xc = M.xc;  yc = M.yc;  hm = D.hom.mesh;
pN = {'\psi_{11}', '\psi_{22}', '\psi_{12}'};
sN = {'\sigma_{11}', '\sigma_{22}', '\sigma_{12}'};
for iN = 1:numel(Ns)
    f = D.het(iN);  hf = f.mesh;  m = M.N(iN);  N = Ns(iN);
    dirN = fullfile(outDir, sprintf('N%d', N));
    for l = 1:numel(loads)
        % displacement: full-scale | micromorphic | homogeneous
        U5 = m.raw.U(:, l);
        ux = {hf.U(1:2:end, l), U5(1:5:end), hm.U(1:2:end, l)};
        uy = {hf.U(2:2:end, l), U5(2:5:end), hm.U(2:2:end, l)};
        nd = {hf.nodes, M.mesh.nodes, hm.nodes};  el = {hf.elems, M.mesh.elems, hm.elems};
        tt = {'full-scale', 'micromorphic', 'homogeneous'};
        limx = vlim(cat(1, ux{:}));  limy = vlim(cat(1, uy{:}));
        fig = new_fig([50 50 1700 620]);
        for k = 1:3
            subplot(2, 3, k);      mesh_panel(gca, nd{k}, el{k}, ux{k}, 'node', ['u_x  ' tt{k}], limx);
            subplot(2, 3, 3 + k);  mesh_panel(gca, nd{k}, el{k}, uy{k}, 'node', ['u_y  ' tt{k}], limy);
        end
        add_title(fig, sprintf('Displacement: N = %d, load: %s', N, loads{l}));
        save_png(fig, fullfile(dirN, sprintf('%s_displacement.png', loads{l})), dpi);

        % psi and sigma on the micromorphic bins: full-scale (smoothed) | micromorphic | difference
        Pf = eval_bins(f.psiCoef(:, :, :, l), xc, yc, Lx, Ly);  Sf = eval_bins(f.sigCoef(:, :, :, l), xc, yc, Lx, Ly);
        Pm = m.Pbin(:, :, :, l);  Sm = m.Sbin(:, :, :, l);
        fig = new_fig([50 50 1900 760]);
        for c = 1:3
            a = Pf(:, :, c);  b = Pm(:, :, c);  d = b - a;  lp = vlim([a(:); b(:)]);  ld = vlim([-max(abs(d(:))); max(abs(d(:)))]);
            subplot(3, 6, (c-1)*6 + 1);  grid_panel(gca, a, Lx, Ly, [pN{c} '  full-scale (smoothed)'], lp);
            subplot(3, 6, (c-1)*6 + 2);  grid_panel(gca, b, Lx, Ly, [pN{c} '  micromorphic'], lp);
            subplot(3, 6, (c-1)*6 + 3);  grid_panel(gca, d, Lx, Ly, 'difference', ld);
            a = Sf(:, :, c);  b = Sm(:, :, c);  d = b - a;  ls = vlim([a(:); b(:)]);  ld = vlim([-max(abs(d(:))); max(abs(d(:)))]);
            subplot(3, 6, (c-1)*6 + 4);  grid_panel(gca, a, Lx, Ly, [sN{c} '  full-scale (smoothed)'], ls);
            subplot(3, 6, (c-1)*6 + 5);  grid_panel(gca, b, Lx, Ly, [sN{c} '  micromorphic'], ls);
            subplot(3, 6, (c-1)*6 + 6);  grid_panel(gca, d, Lx, Ly, 'difference', ld);
        end
        add_title(fig, sprintf('psi (micro strain) and stress, full-scale vs micromorphic: N = %d, load: %s', N, loads{l}));
        save_png(fig, fullfile(dirN, sprintf('%s_psi_sigma.png', loads{l})), dpi);
    end
end
end


function fig_error_summary(C, outDir, dpi)
nL = numel(C.names);  Ns = C.Ns;  col = lines(nL);
fig = new_fig([50 50 1700 520]);
data = {C.errDisp, C.errHomDisp, 'edge displacement error (%)'; ...
        C.errEn, C.errHomEn, 'energy 2U error (%)'; ...
        C.errSig, C.errSigHom, 'stress field misfit (%)'};
for p = 1:3
    ax = subplot(1, 3, p);  hold(ax, 'on');
    labels = cell(1, 2*nL);
    for l = 1:nL
        plot(ax, Ns, data{p, 1}(:, l), '-o', 'Color', col(l, :), 'LineWidth', 1.5);
        plot(ax, Ns, data{p, 2}(:, l), '--s', 'Color', col(l, :), 'LineWidth', 1.2);
        labels{2*l - 1} = [C.names{l} '  micromorphic'];
        labels{2*l} = [C.names{l} '  homogeneous'];
    end
    grid(ax, 'on');  box(ax, 'on');  xlabel(ax, 'N');  ylabel(ax, data{p, 3});
    set(ax, 'XTick', Ns);
    if p == 1, legend(ax, labels, 'Location', 'best', 'FontSize', 8); end
end
add_title(fig, 'Errors with respect to the full-scale plate (test problem)');
save_png(fig, fullfile(outDir, 'error_summary.png'), dpi);
end


%% ========================================================================
%  HELPERS
% ========================================================================

function fig = new_fig(pos)
fig = figure('Visible', 'off', 'Position', pos, 'Color', 'w');
end


function add_title(fig, str)
annotation(fig, 'textbox', [0 0.95 1 0.05], 'String', str, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end


function save_png(fig, file, dpi)
folder = fileparts(file);
if ~isempty(folder) && exist(folder, 'dir') ~= 7, mkdir(folder); end
set(fig, 'Color', 'w');
try
    exportgraphics(fig, file, 'Resolution', dpi);
catch
    print(fig, file, '-dpng', sprintf('-r%d', dpi));
end
close(fig);
end


function mesh_panel(ax, nodes, elems, vals, mode, ttl, lim)
% field on a T6 mesh drawn with its corner triangles; mode 'node': one value per node (smooth), 'elem': one value per element (flat)
vals = double(vals(:));
if strcmp(mode, 'node'), fc = 'interp'; else, fc = 'flat'; end
patch('Parent', ax, 'Faces', double(elems(:, 1:3)), 'Vertices', nodes, 'FaceVertexCData', vals, 'FaceColor', fc, 'EdgeColor', 'none');
axis(ax, 'equal');  axis(ax, 'tight');  box(ax, 'on');
title(ax, ttl, 'FontSize', 9);
if numel(lim) == 2 && all(isfinite(lim)) && lim(2) > lim(1), set(ax, 'CLim', lim); end
colorbar(ax);
end


function grid_panel(ax, Z, Lx, Ly, ttl, lim)
% field on a regular grid of bins (nby x nbx, row 1 = y near 0)
[nby, nbx] = size(Z);
imagesc(ax, [Lx/(2*nbx), Lx*(1 - 1/(2*nbx))], [Ly/(2*nby), Ly*(1 - 1/(2*nby))], double(Z));
set(ax, 'YDir', 'normal');
axis(ax, 'equal');  axis(ax, 'tight');  box(ax, 'on');
title(ax, ttl, 'FontSize', 9);
if numel(lim) == 2 && all(isfinite(lim)) && lim(2) > lim(1), set(ax, 'CLim', lim); end
colorbar(ax);
end


function lim = vlim(v)
% colour limits [min max] of the finite values of v (widened when all values are equal; [] when there are none)
v = double(v(:));  v = v(isfinite(v));
if isempty(v), lim = [];  return, end
lim = [min(v), max(v)];
if lim(2) <= lim(1), lim = lim + [-1 1]*max(abs(lim(1))*1e-6, 1e-12); end
end


function v = elem_mean(Wg, A)
% area-weighted mean over the Gauss points of every element (Wg, A: nElem x nPoints)
Wg = double(Wg);  A = double(A);
v = sum(Wg.*A, 2)./sum(Wg, 2);
end


function Z = eval_bins(C, xc, yc, Lx, Ly)
% smooth Legendre fields C ((degY+1) x (degX+1) x 3) at the bin centres (xc, yc) -> numel(yc) x numel(xc) x 3
Px = legendre_basis(2*xc(:)/Lx - 1, size(C, 2) - 1);  Py = legendre_basis(2*yc(:)/Ly - 1, size(C, 1) - 1);
Z = zeros(numel(yc), numel(xc), 3);
for k = 1:3, Z(:, :, k) = Py*C(:, :, k)*Px.'; end
end
