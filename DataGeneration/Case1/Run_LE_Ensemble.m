clear; clc; close all;

% ============================================================
%  Run_LE_Ensemble.m
%  Drives Main2D_LE.py over a grid of inclusion translations,
%  for a library of linear-elastic load cases (1% strain,
%  single step, single increment).
%
%  Output: ONE ensemble CSV per load case
%          <avgDir>/<LoadCase>_ensemble.csv
%  plus    <avgDir>/Ensemble_Summary.csv  (energy + reactions)
% ============================================================

rootDir = fileparts(mfilename("fullpath"));
addpath(rootDir);

% ------------------------------------------------------------
% helper package
% ------------------------------------------------------------
pkgDir = fullfile(rootDir, "+AllFunctions");
if ~exist(pkgDir,"dir"); mkdir(pkgDir); end

U = "https://raw.githubusercontent.com/mshojaee3/AllFunctionsPub/main/";
F = ["Mat_5A_safeCleanRunDir.m"; ...
     "Mat_7A_updatePyFromParams.m"];
arrayfun(@(f) websave(fullfile(pkgDir,f), U+f), F, 'UniformOutput', false);

% ------------------------------------------------------------
% paths
% ------------------------------------------------------------
srcPy     = fullfile(rootDir, "Main2D_LE.py");
runDir    = fullfile(rootDir, "AbaqusSimulation");
rawDir    = fullfile(rootDir, "LE_Raw_CSV");
energyDir = fullfile(rootDir, "LE_Energy_CSV");

assert(isfile(srcPy), 'Python template not found: %s', srcPy);
for d = [runDir, rawDir, energyDir]
    if ~exist(d,'dir'); mkdir(d); end
end

% ------------------------------------------------------------
% geometry / loading parameters (must be pushed into python)
% ------------------------------------------------------------
geom          = struct();
geom.Lx       = 1.0;
geom.Ly       = 1.0;
geom.Nx       = 1;
geom.Ny       = 1;
geom.Rfrac    = 0.30;
geom.meshFrac = 0.05;

strain0 = 0.01;          % 1 percent

Lx_tot = geom.Nx * geom.Lx;
Ly_tot = geom.Ny * geom.Ly;
R      = geom.Rfrac * geom.Lx;

avgDir = fullfile(rootDir, sprintf('LE_Ensemble_%dx%d', geom.Nx, geom.Ny));
if ~exist(avgDir,'dir'); mkdir(avgDir); end

% ------------------------------------------------------------
% load case library  (names MUST match buildLoadCase in the .py)
% ------------------------------------------------------------
loadCaseList = { ...
    'UniaxialX'; ...
    'UniaxialY'; ...
    'CompressionX'; ...
    'CompressionY'; ...
    'BiaxialTension'; ...
    'BiaxialCompression'; ...
    'BiaxialUnequal'; ...
    'PureShearNormal'; ...
    'SimpleShearXY'; ...
    'SimpleShearYX'; ...
    'PureShearSym'; ...
    'RotationShear'; ...
    'CombinedTensionShear'; ...
    'CombinedBiaxialShear'; ...
    'PureBendingX'; ...
    'PureBendingY'; ...
    'ShearGradientX'; ...
    'ShearGradientY'; ...
    'StretchGradientX'; ...
    'StretchGradientY'; ...
    'DilatationGradient'; ...
    'TensionPlusBending'; ...
    'ShearPlusBending'; ...
    'BiaxialPlusShearGradient'; ...
    'GeneralMixed'; ...
    'FreeUniaxialX'; ...
    'FreeUniaxialY'; ...
    'FreeCompressionY'; ...
    'FreeSimpleShear'; ...
};

% ------------------------------------------------------------
% translation grid
% ------------------------------------------------------------
marginFrac = 0.05;
stepFrac   = 0.05;

zeta1Critical = geom.Lx/2 - R - marginFrac*geom.Lx;
zeta2Critical = geom.Ly/2 - R - marginFrac*geom.Ly;
if zeta1Critical < 0 || zeta2Critical < 0
    error('Critical translation is negative: inclusion or margin too large.');
end

stepX = stepFrac * geom.Lx;
stepY = stepFrac * geom.Ly;
zeta1Vals = -zeta1Critical : stepX : zeta1Critical;
zeta2Vals = -zeta2Critical : stepY : zeta2Critical;

% 1 = diagonal, 2 = full grid, 3 = cross
translationMode = 2;

runList = [];
switch translationMode
    case 1
        nDiag = min(numel(zeta1Vals), numel(zeta2Vals));
        for i = 1:nDiag
            runList = [runList; zeta1Vals(i), zeta2Vals(i)]; %#ok<AGROW>
        end
    case 2
        for i = 1:numel(zeta1Vals)
            for j = 1:numel(zeta2Vals)
                runList = [runList; zeta1Vals(i), zeta2Vals(j)]; %#ok<AGROW>
            end
        end
    case 3
        runList = [0, 0];
        for s = stepX:stepX:zeta1Critical
            runList = [runList; s, 0; -s, 0]; %#ok<AGROW>
        end
        for s = stepY:stepY:zeta2Critical
            runList = [runList; 0, s; 0, -s]; %#ok<AGROW>
        end
    otherwise
        error('Unknown translationMode.');
end

fprintf('\nzeta1Critical = %.4f, zeta2Critical = %.4f\n', zeta1Critical, zeta2Critical);
fprintf('Translations per load case = %d\n', size(runList,1));
fprintf('Load cases = %d, total jobs = %d\n', numel(loadCaseList), ...
        numel(loadCaseList)*size(runList,1));

% ------------------------------------------------------------
% Abaqus
% ------------------------------------------------------------
abaqusBat = 'C:\SIMULIA\Commands\abaqus.bat';
assert(isfile(abaqusBat), 'abaqusBat does not exist: %s', abaqusBat);

% ------------------------------------------------------------
% ensemble grid
% ------------------------------------------------------------
Ngx = 120;
Ngy = 120;
x_grid = linspace(0, Lx_tot, Ngx);
y_grid = linspace(0, Ly_tot, Ngy);
[Xg, Yg] = meshgrid(x_grid, y_grid);

summaryRows = {};

% ============================================================
% MAIN LOOP
% ============================================================
nCases = numel(loadCaseList);
nRuns  = size(runList,1);
jobCounter = 0;

for c = 1:nCases

    loadCaseName = loadCaseList{c};

    fprintf('\n##################################################\n');
    fprintf('Load case %d / %d : %s\n', c, nCases, loadCaseName);
    fprintf('##################################################\n');

    for k = 1:nRuns

        jobCounter = jobCounter + 1;
        z1 = runList(k,1);
        z2 = runList(k,2);

        fprintf('\n[%d / %d] %s   zeta = (%.4f, %.4f)\n', ...
                jobCounter, nCases*nRuns, loadCaseName, z1, z2);

        AllFunctions.Mat_5A_safeCleanRunDir(runDir);

        params              = struct();
        params.Lx           = geom.Lx;
        params.Ly           = geom.Ly;
        params.NX           = geom.Nx;
        params.NY           = geom.Ny;
        params.Rfrac        = geom.Rfrac;
        params.meshFrac     = geom.meshFrac;
        params.strain0      = strain0;
        params.zeta1        = z1;
        params.zeta2        = z2;
        params.loadCaseName = loadCaseName;
        params.OUTDIR       = strrep(rawDir,    '\', '/');
        params.ENERGYDIR    = strrep(energyDir, '\', '/');

        updatedPy = fullfile(runDir, "Main2D_LE_updated.py");
        AllFunctions.Mat_7A_updatePyFromParams(params, srcPy, updatedPy);

        logFile = fullfile(runDir, sprintf('log_%s_zx_%0.3f_zy_%0.3f.txt', ...
                           loadCaseName, z1, z2));
        logFile = strrep(strrep(logFile, '.', 'p'), '-', 'm');

        cmd = sprintf(['cd /d "%s" && call "%s" cae noGUI="%s" > "%s" 2>&1'], ...
                       runDir, abaqusBat, updatedPy, logFile);
        status = system(cmd);

        if status ~= 0
            error(['Abaqus failed: %s, zeta = (%.4f, %.4f).\nSee log:\n%s'], ...
                   loadCaseName, z1, z2, logFile);
        end
    end

    % --------------------------------------------------------
    % ENSEMBLE AVERAGE -> one CSV for this load case
    % --------------------------------------------------------
    caseRawDir = fullfile(rawDir, loadCaseName);
    if ~exist(caseRawDir,'dir')
        warning('Raw folder missing: %s', caseRawDir); continue;
    end

    fileList = dir(fullfile(caseRawDir, '*_Nodal.csv'));
    Nfiles = numel(fileList);
    if Nfiles == 0
        warning('No nodal CSV found in %s', caseRawDir); continue;
    end
    fprintf('\nAveraging %d translations for %s ...\n', Nfiles, loadCaseName);

    % fields carried through the ensemble average
    fieldNames = {'U1','U2','S11','S22','S33','S12','E11','E22','E33','E12'};
    nF = numel(fieldNames);

    fSum  = repmat({zeros(Ngy, Ngx)}, 1, nF);
    fSum2 = repmat({zeros(Ngy, Ngx)}, 1, nF);

    for k = 1:Nfiles
        fname = fullfile(caseRawDir, fileList(k).name);
        T = readtable(fname);

        bad = isnan(T.X) | isnan(T.Y) | isnan(T.U1) | isnan(T.U2);
        T(bad,:) = [];

        for q = 1:nF
            vq = T.(fieldNames{q});
            Vk = griddata(T.X, T.Y, vq, Xg, Yg, 'linear');
            nanmask = isnan(Vk);
            if any(nanmask(:))
                Vk(nanmask) = griddata(T.X, T.Y, vq, Xg(nanmask), Yg(nanmask), 'nearest');
            end
            fSum{q}  = fSum{q}  + Vk;
            fSum2{q} = fSum2{q} + Vk.^2;
        end
    end

    fBar = cell(1,nF);
    fStd = cell(1,nF);
    for q = 1:nF
        fBar{q} = fSum{q} / Nfiles;
        fStd{q} = sqrt(max(fSum2{q}/Nfiles - fBar{q}.^2, 0));
    end

    U1_bar = fBar{1};  U2_bar = fBar{2};
    sigma1 = fStd{1};  sigma2 = fStd{2};

    dx = x_grid(2) - x_grid(1);
    dy = y_grid(2) - y_grid(1);
    [dU1_dX1, dU1_dX2] = gradient(U1_bar, dx, dy);
    [dU2_dX1, dU2_dX2] = gradient(U2_bar, dx, dy);

    F11 = 1 + dU1_dX1;   F12 = dU1_dX2;
    F21 = dU2_dX1;       F22 = 1 + dU2_dX2;
    J   = F11.*F22 - F12.*F21;

    % small-strain measure derived from the averaged displacement field
    % (Eg* = gradient-derived, distinct from the FE strain Ebar_*)
    Eg11 = dU1_dX1;
    Eg22 = dU2_dX2;
    Eg12 = 0.5*(dU1_dX2 + dU2_dX1);

    % ---- assemble the output table column by column ----
    colNames = {'X','Y','U1_bar','U2_bar','sigma_U1','sigma_U2', ...
                'F11','F12','F21','F22','J','Eg11','Eg12','Eg22'};
    colData  = {Xg, Yg, U1_bar, U2_bar, sigma1, sigma2, ...
                F11, F12, F21, F22, J, Eg11, Eg12, Eg22};

    for q = 3:nF   % S11..E12, skipping U1,U2 already written above
        colNames{end+1} = [fieldNames{q} '_bar'];  %#ok<AGROW>
        colData{end+1}  = fBar{q};                 %#ok<AGROW>
        colNames{end+1} = ['sigma_' fieldNames{q}];%#ok<AGROW>
        colData{end+1}  = fStd{q};                 %#ok<AGROW>
    end

    M = zeros(Ngy*Ngx, numel(colData));
    for q = 1:numel(colData)
        M(:,q) = reshape(colData{q}.', [], 1);   % row-major (iy outer, ix inner)
    end

    Tout = array2table(M, 'VariableNames', colNames);
    Tout = addvars(Tout, repmat({loadCaseName}, Ngy*Ngx, 1), ...
                         repmat(Nfiles, Ngy*Ngx, 1), ...
                   'Before', 1, 'NewVariableNames', {'LoadCase','Nfiles'});

    outname = fullfile(avgDir, sprintf('%s_ensemble.csv', loadCaseName));
    writetable(Tout, outname);
    fprintf('Saved ensemble file: %s\n', outname);

    % --------------------------------------------------------
    % average energy / reactions over translations
    % --------------------------------------------------------
    caseEnergyDir = fullfile(energyDir, loadCaseName);
    sumFiles = dir(fullfile(caseEnergyDir, '*_Summary.csv'));
    if ~isempty(sumFiles)
        sumVars = {'ALLSE','RF_LEFT_x','RF_LEFT_y','RF_RIGHT_x','RF_RIGHT_y', ...
                   'RF_BOTTOM_x','RF_BOTTOM_y','RF_TOP_x','RF_TOP_y','Volume', ...
                   'Savg_11','Savg_22','Savg_33','Savg_12', ...
                   'Eavg_11','Eavg_22','Eavg_33','Eavg_12','Uavg_1','Uavg_2'};
        Sall = [];
        for k = 1:numel(sumFiles)
            Tk = readtable(fullfile(caseEnergyDir, sumFiles(k).name));
            Sall = [Sall; Tk{1, sumVars}]; %#ok<AGROW>
        end
        row = [{loadCaseName, numel(sumFiles), mean(Sall(:,1)), std(Sall(:,1))}, ...
               num2cell(mean(Sall(:,2:end), 1))];
        summaryRows(end+1,:) = row; %#ok<AGROW>
    end
end

% ============================================================
% GLOBAL SUMMARY
% ============================================================
if ~isempty(summaryRows)
    Tsum = cell2table(summaryRows, 'VariableNames', ...
        {'LoadCase','Ntrans','ALLSE_mean','ALLSE_std', ...
         'RF_LEFT_x','RF_LEFT_y','RF_RIGHT_x','RF_RIGHT_y', ...
         'RF_BOTTOM_x','RF_BOTTOM_y','RF_TOP_x','RF_TOP_y','Volume', ...
         'Savg_11','Savg_22','Savg_33','Savg_12', ...
         'Eavg_11','Eavg_22','Eavg_33','Eavg_12','Uavg_1','Uavg_2'});
    writetable(Tsum, fullfile(avgDir, 'Ensemble_Summary.csv'));
    fprintf('\nSaved global summary: %s\n', fullfile(avgDir,'Ensemble_Summary.csv'));
end

fprintf('\nAll load cases finished.\n');

% ------------------------------------------------------------
% quick check of the last load case
% ------------------------------------------------------------
figure('Name', ['Ensemble F: ' loadCaseName], 'Position',[100 100 1000 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
comps = {F11, F12, F21, F22};
names = {'F_{11}','F_{12}','F_{21}','F_{22}'};
for q = 1:4
    nexttile;
    imagesc(x_grid, y_grid, comps{q}); axis xy equal tight; colorbar;
    xlabel('X_1'); ylabel('X_2'); title(names{q});
end
sgtitle(sprintf('Ensemble average, load case: %s', loadCaseName), 'Interpreter','none');
