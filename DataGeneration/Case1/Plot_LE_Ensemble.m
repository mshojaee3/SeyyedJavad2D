clear; clc; close all;

% ============================================================
%  Plot_LE_Ensemble.m
%  Standalone viewer for the ensemble CSVs written by
%  Run_LE_Ensemble.m:
%       <avgDir>/<LoadCase>_ensemble.csv
%
%  Reads one load case, rebuilds the structured grid from the
%  X/Y columns, and plots displacement, ensemble fluctuation,
%  deformation gradient, FE stress/strain and line profiles.
% ============================================================

% ------------------------------------------------------------
% USER INPUT
% ------------------------------------------------------------
loadCaseName = 'UniaxialX';      % which case to plot
Nx = 1;  Ny = 1;                 % must match the run (folder suffix)

savePNG   = false;               % write PNGs next to the CSV
dispScale = 20;                  % magnification of the deformed overlay
nQuiver   = 20;                  % arrows per direction

% ------------------------------------------------------------
% locate the file
% ------------------------------------------------------------
rootDir = fileparts(mfilename("fullpath"));
avgDir  = fullfile(rootDir, sprintf('LE_Ensemble_%dx%d', Nx, Ny));
csvFile = fullfile(avgDir, sprintf('%s_ensemble.csv', loadCaseName));

if ~isfile(csvFile)
    d = dir(fullfile(avgDir, '*_ensemble.csv'));
    if isempty(d)
        error('No ensemble CSV found in:\n%s', avgDir);
    end
    fprintf('Available load cases in %s:\n', avgDir);
    for k = 1:numel(d)
        fprintf('   %s\n', erase(d(k).name, '_ensemble.csv'));
    end
    error('File not found: %s', csvFile);
end

T = readtable(csvFile);
fprintf('Loaded %s\n  %d rows, %d translations averaged\n', ...
        csvFile, height(T), T.Nfiles(1));

% ------------------------------------------------------------
% rebuild the structured grid (independent of row ordering)
% ------------------------------------------------------------
[xu, ~, ixv] = unique(T.X);
[yu, ~, iyv] = unique(T.Y);
Ngx = numel(xu);
Ngy = numel(yu);

x_grid = xu(:).';
y_grid = yu(:).';
[Xg, Yg] = meshgrid(x_grid, y_grid);

% anonymous reshaper: scattered rows -> (Ngy x Ngx) matrix
G = @(name) accumarray([iyv, ixv], T.(name), [Ngy, Ngx], @mean, NaN);

fprintf('  grid %d x %d, domain [%g, %g] x [%g, %g]\n', ...
        Ngx, Ngy, min(xu), max(xu), min(yu), max(yu));

% ------------------------------------------------------------
% pull fields
% ------------------------------------------------------------
U1  = G('U1_bar');     U2  = G('U2_bar');
sU1 = G('sigma_U1');   sU2 = G('sigma_U2');

F11 = G('F11');   F12 = G('F12');
F21 = G('F21');   F22 = G('F22');
Jdet = G('J');

Eg11 = G('Eg11');  Eg12 = G('Eg12');  Eg22 = G('Eg22');

S11 = G('S11_bar');  S22 = G('S22_bar');
S33 = G('S33_bar');  S12 = G('S12_bar');
sS11 = G('sigma_S11');

% Abaqus reports engineering shear in E12, so eps12 = E12/2
E11 = G('E11_bar');  E22 = G('E22_bar');
E33 = G('E33_bar');  E12 = 0.5*G('E12_bar');
sE11 = G('sigma_E11');

Umag = sqrt(U1.^2 + U2.^2);
vM   = sqrt(0.5*((S11-S22).^2 + (S22-S33).^2 + (S33-S11).^2 + 6*S12.^2));

% ------------------------------------------------------------
% console summary
% ------------------------------------------------------------
mn = @(A) mean(A(:), 'omitnan');
fprintf('\n--- domain means, %s ---\n', loadCaseName);
fprintf('  <u1>  = %+.6e     <u2>  = %+.6e\n', mn(U1), mn(U2));
fprintf('  <F>   = [%+.6f %+.6f ; %+.6f %+.6f],   <J> = %.6f\n', ...
        mn(F11), mn(F12), mn(F21), mn(F22), mn(Jdet));
fprintf('  <s11> = %+.6e  <s22> = %+.6e  <s33> = %+.6e  <s12> = %+.6e\n', ...
        mn(S11), mn(S22), mn(S33), mn(S12));
fprintf('  <e11> = %+.6e  <e22> = %+.6e  <e33> = %+.6e  <e12> = %+.6e\n', ...
        mn(E11), mn(E22), mn(E33), mn(E12));
fprintf('  FE vs gradient:  <e11>-<Eg11> = %+.3e    <e12>-<Eg12> = %+.3e\n', ...
        mn(E11)-mn(Eg11), mn(E12)-mn(Eg12));
fprintf('  max fluctuation: sigma_u1 = %.3e   sigma_u2 = %.3e\n', ...
        max(sU1(:)), max(sU2(:)));

figs = gobjects(0);

% ============================================================
% FIG 1 - displacement
% ============================================================
figs(end+1) = figure('Name',[loadCaseName '_displacement'], ...
                     'Position',[60 60 1100 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
panel(x_grid, y_grid, U1,   'mean u_1');
panel(x_grid, y_grid, U2,   'mean u_2');
panel(x_grid, y_grid, Umag, '|mean u|');

nexttile;
sx = max(1, round(Ngx/nQuiver));
sy = max(1, round(Ngy/nQuiver));
quiver(Xg(1:sy:end,1:sx:end), Yg(1:sy:end,1:sx:end), ...
       U1(1:sy:end,1:sx:end), U2(1:sy:end,1:sx:end), 0, 'k'); hold on;
plot(Xg(1:sy:end,1:sx:end) + dispScale*U1(1:sy:end,1:sx:end), ...
     Yg(1:sy:end,1:sx:end) + dispScale*U2(1:sy:end,1:sx:end), '.', ...
     'Color',[0.85 0.33 0.10], 'MarkerSize', 6);
axis equal tight; grid on; xlabel('X_1'); ylabel('X_2');
title(sprintf('displacement (deformed markers x%g)', dispScale));

% ============================================================
% FIG 2 - ensemble fluctuation across translations
% ============================================================
figs(end+1) = figure('Name',[loadCaseName '_fluctuation'], ...
                     'Position',[80 80 1100 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
panel(x_grid, y_grid, sU1,  'std of u_1');
panel(x_grid, y_grid, sU2,  'std of u_2');
panel(x_grid, y_grid, sS11, 'std of \sigma_{11}');
panel(x_grid, y_grid, sE11, 'std of \epsilon_{11}');

% ============================================================
% FIG 3 - deformation gradient
% ============================================================
figs(end+1) = figure('Name',[loadCaseName '_F'], ...
                     'Position',[100 100 1100 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
panel(x_grid, y_grid, F11, 'F_{11}');
panel(x_grid, y_grid, F12, 'F_{12}');
panel(x_grid, y_grid, F21, 'F_{21}');
panel(x_grid, y_grid, F22, 'F_{22}');

% ============================================================
% FIG 4 - stress
% ============================================================
figs(end+1) = figure('Name',[loadCaseName '_stress'], ...
                     'Position',[120 120 1100 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
panel(x_grid, y_grid, S11, '\sigma_{11}');
panel(x_grid, y_grid, S22, '\sigma_{22}');
panel(x_grid, y_grid, S12, '\sigma_{12}');
panel(x_grid, y_grid, vM,  'von Mises');

% ============================================================
% FIG 5 - FE strain vs gradient-derived strain
% ============================================================
figs(end+1) = figure('Name',[loadCaseName '_strain'], ...
                     'Position',[140 140 1300 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
panel(x_grid, y_grid, E11,  '\epsilon_{11} (FE nodal)');
panel(x_grid, y_grid, E22,  '\epsilon_{22} (FE nodal)');
panel(x_grid, y_grid, E12,  '\epsilon_{12} (FE nodal)');
panel(x_grid, y_grid, Eg11, '\epsilon_{11} from grad(u)');
panel(x_grid, y_grid, Eg22, '\epsilon_{22} from grad(u)');
panel(x_grid, y_grid, Eg12, '\epsilon_{12} from grad(u)');

% ============================================================
% FIG 6 - line profiles through the centre
% ============================================================
iyMid = max(1, round(Ngy/2));
ixMid = max(1, round(Ngx/2));

figs(end+1) = figure('Name',[loadCaseName '_profiles'], ...
                     'Position',[160 160 1200 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(x_grid, U1(iyMid,:), 'LineWidth',1.4); hold on;
plot(x_grid, U2(iyMid,:), 'LineWidth',1.4);
grid on; xlabel('X_1'); ylabel('displacement');
legend('u_1','u_2','Location','best');
title(sprintf('horizontal cut at X_2 = %.3f', y_grid(iyMid)));

nexttile;
plot(y_grid, U1(:,ixMid), 'LineWidth',1.4); hold on;
plot(y_grid, U2(:,ixMid), 'LineWidth',1.4);
grid on; xlabel('X_2'); ylabel('displacement');
legend('u_1','u_2','Location','best');
title(sprintf('vertical cut at X_1 = %.3f', x_grid(ixMid)));

nexttile;
plot(x_grid, S11(iyMid,:), 'LineWidth',1.4); hold on;
plot(x_grid, S22(iyMid,:), 'LineWidth',1.4);
plot(x_grid, S12(iyMid,:), 'LineWidth',1.4);
grid on; xlabel('X_1'); ylabel('stress');
legend('\sigma_{11}','\sigma_{22}','\sigma_{12}','Location','best');
title('stress along horizontal cut');

nexttile;
plot(x_grid, E11(iyMid,:), 'LineWidth',1.4); hold on;
plot(x_grid, Eg11(iyMid,:), '--', 'LineWidth',1.4);
grid on; xlabel('X_1'); ylabel('\epsilon_{11}');
legend('FE nodal','from grad(u)','Location','best');
title('FE strain vs gradient of averaged displacement');

% ============================================================
% super-titles and optional export
% ============================================================
for k = 1:numel(figs)
    figure(figs(k));
    sgtitle(sprintf('%s   -   ensemble of %d translations', ...
            loadCaseName, T.Nfiles(1)), 'Interpreter','none');
end

if savePNG
    for k = 1:numel(figs)
        exportgraphics(figs(k), ...
            fullfile(avgDir, [figs(k).Name '.png']), 'Resolution', 200);
    end
    fprintf('\nPNG files written to %s\n', avgDir);
end

% ============================================================
% LOCAL FUNCTIONS (must stay at the end of the script)
% ============================================================
function panel(xg, yg, A, ttl)
    nexttile;
    imagesc(xg, yg, A);
    axis xy equal tight;
    colorbar;
    xlabel('X_1'); ylabel('X_2');
    title(ttl);
end