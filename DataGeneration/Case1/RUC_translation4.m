clear; clc; close all;

% Root directory = folder of this script
rootDir = fileparts(mfilename("fullpath"));
addpath(rootDir);

% Packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkgDir  = fullfile(rootDir, "+AllFunctions");
if ~exist(pkgDir,"dir"); mkdir(pkgDir); end

U = "https://raw.githubusercontent.com/mshojaee3/AllFunctionsPub/main/";
F = ["Mat_5A_safeCleanRunDir.m"; ...
     "Mat_7A_updatePyFromParams.m"];

arrayfun(@(f) websave(fullfile(pkgDir,f), U+f), F, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------
% paths
% ------------------------------------------------------------
srcPy  = fullfile(rootDir, "Main2D_II.py");
runDir = fullfile(rootDir, "AbaqusSimulation");
outDir = fullfile(rootDir, "Translated_RUC_AllFrames_CSV");

if ~exist(runDir,'dir')
    mkdir(runDir);
end
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% ------------------------------------------------------------
% geometry parameters
% IMPORTANT: these must match the values inside Main2D_II.py
% ------------------------------------------------------------
geom = struct();
geom.Lx = 1.0;
geom.Ly = 1.0;
geom.Nx = 6;
geom.Ny = 6;

% ------------------------------------------------------------
% output folders that depend on RUC number
% ------------------------------------------------------------
energyFolderName = sprintf('Translated_RUC_Energy_%dx%d', geom.Nx, geom.Ny);
energyDir = fullfile(rootDir, energyFolderName);
if ~exist(energyDir,'dir')
    mkdir(energyDir);
end

avgFolderName = sprintf('Translated_RUC_Ensemble_%dx%d', geom.Nx, geom.Ny);
avgDir = fullfile(rootDir, avgFolderName);
if ~exist(avgDir,'dir')
    mkdir(avgDir);
end

% Inclusion radius
R = 0.3 * geom.Lx;

% ------------------------------------------------------------
% load cases
% ------------------------------------------------------------
loadCaseList = { ...
    'PureBending'; ...
    'TensionBiaxial'; ...
    'CompressiveUniaxial'; ...
    'CompressiveBiaxial'; ...
    'SimpleShear'; ...
    'PureShear'; ...
    'ParabolicBending'; ...
    'TensionUniaxial'; ...
};

% ------------------------------------------------------------
% translation settings
% marginFrac = required remaining gap between inclusion edge
% and RUC boundary, measured as fraction of RUC length
% ------------------------------------------------------------
marginFrac = 0.05;   % marginFrac = 0.05  --> remaining gap = 5% of L
stepFrac   = 0.05;   % translation step = 5% of each RUC length

% critical translations
zeta1Critical = geom.Lx/2 - R - marginFrac*geom.Lx;
zeta2Critical = geom.Ly/2 - R - marginFrac*geom.Ly;

if zeta1Critical < 0 || zeta2Critical < 0
    error(['Critical translation is negative. ' ...
           'The inclusion is too large or the chosen margin is too large.']);
end

stepX = stepFrac * geom.Lx;
stepY = stepFrac * geom.Ly;

zeta1Vals = -zeta1Critical : stepX : zeta1Critical;
zeta2Vals = -zeta2Critical : stepY : zeta2Critical;

% ------------------------------------------------------------
% translation mode
% 1 = diagonal only
% 2 = full grid
% 3 = cross only (center + right/left/up/down)
% ------------------------------------------------------------
translationMode = 2;

% ------------------------------------------------------------
% build run list
% ------------------------------------------------------------
runList = [];

if translationMode == 1
    nDiag = min(numel(zeta1Vals), numel(zeta2Vals));
    for i = 1:nDiag
        runList = [runList; zeta1Vals(i), zeta2Vals(i)];
    end

elseif translationMode == 2
    for i = 1:numel(zeta1Vals)
        for j = 1:numel(zeta2Vals)
            runList = [runList; zeta1Vals(i), zeta2Vals(j)];
        end
    end

elseif translationMode == 3
    runList = [0, 0];

    % right
    for s = stepX:stepX:zeta1Critical
        runList = [runList; s, 0];
    end

    % left
    for s = stepX:stepX:zeta1Critical
        runList = [runList; -s, 0];
    end

    % up
    for s = stepY:stepY:zeta2Critical
        runList = [runList; 0, s];
    end

    % down
    for s = stepY:stepY:zeta2Critical
        runList = [runList; 0, -s];
    end

else
    error('Unknown translationMode.');
end

fprintf('\nCritical translations:\n');
fprintf('zeta1Critical = %.4f\n', zeta1Critical);
fprintf('zeta2Critical = %.4f\n', zeta2Critical);
fprintf('Total runs      = %d\n', size(runList,1));

% ------------------------------------------------------------
% Abaqus path
% ------------------------------------------------------------
abaqusBat = 'C:\SIMULIA\Commands\abaqus.bat';
assert(isfile(abaqusBat), 'abaqusBat does not exist: %s', abaqusBat);

% ------------------------------------------------------------
% ensemble averaging settings
% ------------------------------------------------------------
Ngx = 120;
Ngy = 120;

Lx_tot = geom.Nx * geom.Lx;
Ly_tot = geom.Ny * geom.Ly;

x_grid = linspace(0, Lx_tot, Ngx);
y_grid = linspace(0, Ly_tot, Ngy);
[Xg, Yg] = meshgrid(x_grid, y_grid);

% ------------------------------------------------------------
% loop over translations
% averaging is done immediately after each load case finishes
% ------------------------------------------------------------
nCases = numel(loadCaseList);
nRuns  = size(runList,1);
totalJobs = nCases * nRuns;
jobCounter = 0;

for c = 1:nCases

    loadCaseName = loadCaseList{c};

    fprintf('\n##################################################\n');
    fprintf('Starting load case %d / %d : %s\n', c, nCases, loadCaseName);
    fprintf('##################################################\n');

    for k = 1:nRuns

        jobCounter = jobCounter + 1;

        z1 = runList(k,1);
        z2 = runList(k,2);

        fprintf('\n========================================\n');
        fprintf('Global job %d / %d\n', jobCounter, totalJobs);
        fprintf('Load case   : %s\n', loadCaseName);
        fprintf('Run %d / %d\n', k, nRuns);
        fprintf('zeta1 = %.4f, zeta2 = %.4f\n', z1, z2);
        fprintf('========================================\n');

        % clean working directory
        AllFunctions.Mat_5A_safeCleanRunDir(runDir);

        % parameters to overwrite in python template
        params = struct();
        params.Lx           = geom.Lx;
        params.Ly           = geom.Ly;
        params.NX           = geom.Nx;
        params.NY           = geom.Ny;
        params.zeta1        = z1;
        params.zeta2        = z2;
        params.loadCaseName = loadCaseName;
        params.OUTDIR       = strrep(outDir, '\', '/');
        params.ENERGYDIR    = strrep(energyDir, '\', '/');

        updatedPy = fullfile(runDir, "Main2D_II_updated.py");

        % create updated python file
        AllFunctions.Mat_7A_updatePyFromParams(params, srcPy, updatedPy);

        % log file name per load case and translation
        logFile = fullfile(runDir, sprintf('log_%s_zx_%0.3f_zy_%0.3f.txt', ...
            loadCaseName, z1, z2));
        logFile = strrep(logFile, '.', 'p');
        logFile = strrep(logFile, '-', 'm');

        % run Abaqus directly
        cmd = sprintf(['cd /d "%s" && ' ...
                       'call "%s" cae noGUI="%s" > "%s" 2>&1'], ...
                       runDir, abaqusBat, updatedPy, logFile);

        disp("CMD = " + string(cmd));

        status = system(cmd);

        if status ~= 0
            error(['Abaqus failed for load case %s, zeta1=%.4f, zeta2=%.4f.\n' ...
                   'Check log file:\n%s'], ...
                   loadCaseName, z1, z2, logFile);
        end
    end

    fprintf('\nFinished all translations for load case: %s\n', loadCaseName);
    fprintf('Starting ensemble averaging immediately for this load case...\n');

    % --------------------------------------------------------
    % average this load case now
    % --------------------------------------------------------
    caseRawDir = fullfile(outDir, loadCaseName);
    caseAvgDir = fullfile(avgDir, loadCaseName);

    if ~exist(caseRawDir, 'dir')
        warning('Raw folder does not exist: %s', caseRawDir);
        continue;
    end

    if ~exist(caseAvgDir, 'dir')
        mkdir(caseAvgDir);
    end

    % find all raw csv files of this load case
    fileList = dir(fullfile(caseRawDir, '*.csv'));

    if isempty(fileList)
        warning('No CSV files found in %s', caseRawDir);
        continue;
    end

    % --------------------------------------------------------
    % detect all frame IDs from filenames
    % expects names like ..._f0000_Nodal.csv
    % --------------------------------------------------------
    frameIDs = [];

    for k = 1:numel(fileList)
        fname = fileList(k).name;
        tok = regexp(fname, '_f(\d+)_Nodal\.csv$', 'tokens');
        if ~isempty(tok)
            frameIDs(end+1) = str2double(tok{1}{1});
        end
    end

    frameIDs = unique(frameIDs);

    if isempty(frameIDs)
        warning('No frame IDs detected in %s', caseRawDir);
        continue;
    end

    fprintf('Found %d frame IDs for %s\n', numel(frameIDs), loadCaseName);

    for iFrame = 1:numel(frameIDs)

        frameID = frameIDs(iFrame);

        fprintf('\nProcessing %s frame f%04d\n', loadCaseName, frameID);

        % pattern to collect all translations of this frame
        pattern = sprintf('*_f%04d_Nodal.csv', frameID);
        frameFiles = dir(fullfile(caseRawDir, pattern));
        Nfiles = numel(frameFiles);

        if Nfiles == 0
            warning('No files found for %s frame f%04d', loadCaseName, frameID);
            continue;
        end

        fprintf('Found %d translation files for this frame.\n', Nfiles);

        % allocate
        U1_sum  = zeros(Ngy, Ngx);
        U2_sum  = zeros(Ngy, Ngx);
        U1_sum2 = zeros(Ngy, Ngx);
        U2_sum2 = zeros(Ngy, Ngx);

        stepTime_ref = NaN;

        % ----------------------------------------------------
        % loop over all translations of this frame
        % ----------------------------------------------------
        for k = 1:Nfiles

            fname = fullfile(caseRawDir, frameFiles(k).name);
            fprintf('  Reading %d / %d : %s\n', k, Nfiles, frameFiles(k).name);

            T = readtable(fname);

            if ismember('StepTime', T.Properties.VariableNames)
                stepTime_ref = T.StepTime(1);
            end

            bad = isnan(T.X) | isnan(T.Y) | isnan(T.U1) | isnan(T.U2);
            T(bad,:) = [];

            x  = T.X;
            y  = T.Y;
            u1 = T.U1;
            u2 = T.U2;

            U1k = griddata(x, y, u1, Xg, Yg, 'linear');
            U2k = griddata(x, y, u2, Xg, Yg, 'linear');

            nanmask = isnan(U1k) | isnan(U2k);
            if any(nanmask(:))
                U1k(nanmask) = griddata(x, y, u1, Xg(nanmask), Yg(nanmask), 'nearest');
                U2k(nanmask) = griddata(x, y, u2, Xg(nanmask), Yg(nanmask), 'nearest');
            end

            U1_sum  = U1_sum  + U1k;
            U2_sum  = U2_sum  + U2k;
            U1_sum2 = U1_sum2 + U1k.^2;
            U2_sum2 = U2_sum2 + U2k.^2;
        end

        % ----------------------------------------------------
        % ensemble mean and fluctuation
        % ----------------------------------------------------
        U1_bar = U1_sum / Nfiles;
        U2_bar = U2_sum / Nfiles;

        U1_var = U1_sum2 / Nfiles - U1_bar.^2;
        U2_var = U2_sum2 / Nfiles - U2_bar.^2;

        U1_var = max(U1_var, 0);
        U2_var = max(U2_var, 0);

        sigma1 = sqrt(U1_var);
        sigma2 = sqrt(U2_var);

        % ----------------------------------------------------
        % deformation gradient from ensemble mean
        % ----------------------------------------------------
        dx = x_grid(2) - x_grid(1);
        dy = y_grid(2) - y_grid(1);

        [dU1_dX1, dU1_dX2] = gradient(U1_bar, dx, dy);
        [dU2_dX1, dU2_dX2] = gradient(U2_bar, dx, dy);

        F11 = 1 + dU1_dX1;
        F12 =     dU1_dX2;
        F21 =     dU2_dX1;
        F22 = 1 + dU2_dX2;

        J = F11 .* F22 - F12 .* F21;

        % ----------------------------------------------------
        % save ensemble csv
        % ----------------------------------------------------
        outname = fullfile(caseAvgDir, sprintf('%s_f%04d_u_bar.csv', loadCaseName, frameID));

        fid = fopen(outname, 'w');
        fprintf(fid, 'FrameID,StepTime,X,Y,U1_bar,U2_bar,sigma1,sigma2,F11,F12,F21,F22,J\n');

        for iy = 1:Ngy
            for ix = 1:Ngx
                fprintf(fid, '%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
                    frameID, stepTime_ref, ...
                    Xg(iy,ix), Yg(iy,ix), ...
                    U1_bar(iy,ix), U2_bar(iy,ix), ...
                    sigma1(iy,ix), sigma2(iy,ix), ...
                    F11(iy,ix), F12(iy,ix), F21(iy,ix), F22(iy,ix), J(iy,ix));
            end
        end
        fclose(fid);

        fprintf('Saved ensemble file: %s\n', outname);
    end

    fprintf('Finished ensemble averaging for load case: %s\n', loadCaseName);
end

fprintf('\nAll load cases, translations, and ensemble averaging finished.\n');




%  optional figure for F components of the last saved data
figure('Name','Last averaged F components', ...
       'Position',[100 100 1000 800]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(x_grid, y_grid, F11); axis xy equal tight; colorbar;
xlabel('X_1'); ylabel('X_2');
title('F_{11}');

nexttile;
imagesc(x_grid, y_grid, F12); axis xy equal tight; colorbar;
xlabel('X_1'); ylabel('X_2');
title('F_{12}');

nexttile;
imagesc(x_grid, y_grid, F21); axis xy equal tight; colorbar;
xlabel('X_1'); ylabel('X_2');
title('F_{21}');

nexttile;
imagesc(x_grid, y_grid, F22); axis xy equal tight; colorbar;
xlabel('X_1'); ylabel('X_2');
title('F_{22}');

sgtitle(sprintf('Last averaged result: %s, frame f%04d', loadCaseName, frameID));