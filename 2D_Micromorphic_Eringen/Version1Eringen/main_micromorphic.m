function main_micromorphic(quick, caseList, force)
% =========================================================================
%  MAIN WORKFLOW
%
%  CASE 1 = calibration / identification
%  CASE 2 = validation only
%
%  Workflow:
%
%                    CASE 1
%                       |
%                       v
%             full-scale identification
%                       |
%                       v
%                 identification
%                       |
%                       v
%                      par1
%                       |
%                 +-----+-----+
%                 |           |
%                 v           v
%             CASE 2       CASE 2
%          full-scale    micromorphic
%             test          test
%                 |           |
%                 +-----+-----+
%                       |
%                       v
%                  comparison
%
%  IMPORTANT:
%  - C_hom is calculated once from homogenize2D_PBC.
%  - The same C_hom is used in CASE 1 and CASE 2.
%  - Micromorphic parameters are identified ONCE from CASE 1.
%  - The identified parameter set par1 is used unchanged in CASE 2.
%
%  Usage:
%     main_micromorphic                              all stages, reusing everything that is unchanged
%     main_micromorphic(true)                        quick check (own folder results_quick)
%     main_micromorphic(false, {'case2'})
%     main_micromorphic(false, {'case2'}, {'identification'})   recompute the identification (and what depends on it)
%     main_micromorphic(false, {'case2'}, 'all')                recompute everything
%
%  SAVING AND REUSE
%  Every stage saves its result together with the input it was made with.  When you run again, a stage whose input (parameters.m values) is
%  unchanged is LOADED instead of recomputed; if you change a value, only the stages that depend on it are recomputed.
%  Stages that can be forced (third argument, a name or a cell array of names; the stages that depend on them are forced too):
%     'homog'  'fullIdent'  'identification'  'fullTest'  'micro'      ('all' = everything)
%  After changing the CODE (not parameters.m), increase P.out.version in parameters.m, or force the stages concerned.
%  The comparison is cheap and always recomputed from the saved data.
%
%  FOLDERS  (root = results, or results_quick in a quick run)
%     homogenization/                  homogenization.mat (C_hom, Cs33, displacement + Gauss-point strain/stress of the 3 unit-strain cases),
%                                      C_hom.csv, log.txt, figures/
%     case1/fullscale/                 hom_plate.mat, het_N<N>.mat   full-scale data of the identification problem: nodal displacement, Gauss-point
%                                      strain/stress, bin means before smoothing, Legendre-smoothed fields, coefficients
%     case1/identification/            identification.mat, par1.mat
%     case1/calibration.mat            P1, Dident (without the large mesh fields), par1, C_hom, Cs33
%     case2/fullscale/                 same as case1/fullscale for the test problem
%     case2/micromorphic/              micromorphic.mat   micromorphic solution: nodal unknowns, Gauss-point strain/psi/stress, bin means
%     case2/results.mat, comparison.csv
%     caseN/logs/                      screen output of every computed stage
%     */figures/                       PNG files, written at the end of the run (P.out.png)
% =========================================================================

if nargin < 1 || isempty(quick), quick = false; end
if nargin < 2 || isempty(caseList), caseList = {'case2'}; end
if nargin < 3 || isempty(force), force = {}; end
if ischar(force) || isstring(force), force = cellstr(force); end
fz = expand_force(force);
tAll = tic;


%% ========================================================================
%  GLOBAL HOMOGENIZED CLASSICAL STIFFNESS
% ========================================================================

P1 = parameters('case1', quick);
reuse = P1.out.reuse;
homDir = P1.out.homFolder;
dir1 = P1.out.folder;

fprintf('\n');
fprintf('############################################################\n');
fprintf('# HOMOGENIZED CLASSICAL STIFFNESS C_hom                    #\n');
fprintf('############################################################\n');

keyH = struct('version', P1.out.version, 'material', rmfield(P1.material, 'Chom'), 'meshFrac', P1.fe.meshFrac);
H = stage_cache(fullfile(homDir, 'homogenization.mat'), keyH, reuse && ~fz.homog, ...
                @() run_homog(quick), fullfile(homDir, 'log.txt'));
C_hom = H.C_hom;
Cs33 = H.Cs33;
P1.material.Chom = C_hom;
write_matrix(C_hom, fullfile(homDir, 'C_hom.csv'));

fprintf('\nC_hom used by the micromorphic model:\n');
disp(C_hom);
fprintf('Cs33:\n');
disp(Cs33);


%% ========================================================================
%  CASE 1: CALIBRATION / IDENTIFICATION
% ========================================================================

fprintf('\n');
fprintf('################################################################\n');
fprintf('################  CASE 1  CALIBRATION  ########################\n');
fprintf('################################################################\n');

%% STEP 1: full-scale identification data

fprintf('\nSTEP 1  full-scale data of the identification problem\n');

Dident = fullscale_solution(P1, 'ident', fullfile(dir1, 'fullscale'), fz.fullIdent);


%% STEP 2: identify micromorphic parameters

fprintf('\n');
fprintf('STEP 2  identification of the micromorphic parameters\n');

keyI = struct( ...
    'version', P1.out.version, ...
    'material', P1.material, ...
    'ident', P1.ident, ...
    'model', P1.model, ...
    'residual', P1.residual, ...
    'fe', struct('meshFrac', P1.fe.meshFrac, 'NB', P1.fe.NB, 'psiDeg', P1.fe.psiDeg, 'nyCalib', P1.fe.nyCalib));

par1 = stage_cache(fullfile(dir1, 'identification', 'identification.mat'), keyI, reuse && ~fz.identification, ...
                   @() identification(P1, Dident), fullfile(dir1, 'logs', 'identification.txt'));
save(fullfile(dir1, 'identification', 'par1.mat'), 'par1');

fprintf('\n');
fprintf('CASE 1 complete: parameter set par1 identified.\n');


%% Save Case 1 calibration (light version; the large mesh fields are in case1/fullscale)

s = struct('P1', P1, 'Dident', strip_fields(Dident), 'par1', par1, 'C_hom', C_hom, 'Cs33', Cs33);
save(fullfile(dir1, 'calibration.mat'), '-struct', 's');


%% ========================================================================
%  CASE 2: VALIDATION ONLY
% ========================================================================

res = cell(1, numel(caseList));

for ic = 1:numel(caseList)

    caseName = caseList{ic};

    if ~strcmpi(caseName, 'case2')
        error( ...
            'Validation case ''%s'' is not supported. Use case2.', ...
            caseName);
    end


    %% Load Case 2 parameters

    P2 = parameters('case2', quick);

    % Use EXACTLY the same homogenized classical stiffness calculated
    % from the material definition above.
    P2.material.Chom = C_hom;

    dir2 = P2.out.folder;


    fprintf('\n');
    fprintf('################################################################\n');
    fprintf('################  CASE 2  VALIDATION  #########################\n');
    fprintf('################################################################\n');

    fprintf('\nUsing Case 1 identified parameter set par1.\n');
    fprintf('NO identification is performed for Case 2.\n');

    fprintf('\nC_hom used for Case 2:\n');
    disp(P2.material.Chom);


    %% STEP 3: full-scale Case 2 test problem

    fprintf('\n');
    fprintf('STEP 3  full-scale solution of the Case 2 test problem\n');

    Dtest = fullscale_solution(P2, 'test', fullfile(dir2, 'fullscale'), fz.fullTest);


    %% STEP 4: Case 2 micromorphic prediction
    %
    % IMPORTANT:
    % par1 comes ONLY from Case 1.
    % There is no call to identification(...) here.

    fprintf('\n');
    fprintf('STEP 4  micromorphic solution of the Case 2 test problem\n');

    keyM = struct( ...
        'version', P2.out.version, ...
        'par', par1, ...
        'material', P2.material, ...
        'test', P2.test, ...
        'fe', struct('nyMesh', P2.fe.nyMesh, 'nBinY', P2.fe.nBinY));

    M = stage_cache(fullfile(dir2, 'micromorphic', 'micromorphic.mat'), keyM, reuse && ~fz.micro, ...
                    @() micromorphic_solution(P2, par1), fullfile(dir2, 'logs', 'micromorphic.txt'));


    %% STEP 5: comparison (cheap: always recomputed from the saved data)

    fprintf('\n');
    fprintf('STEP 5  comparison: Case 2 full-scale vs Case 2 micromorphic\n');

    lg = start_log(fullfile(dir2, 'logs', 'comparison.txt'));
    C = compare_results(P2, Dtest, M);
    clear lg

    write_comparison_csv(C, Dtest, fullfile(dir2, 'comparison.csv'));


    %% Save Case 2 results (light version; the large fields are in case2/fullscale and case2/micromorphic)

    Mlite = rmfield(M, 'mesh');
    Mlite.N = rmfield(M.N, 'raw');

    s = struct('P2', P2, 'par1', par1, 'Dtest', strip_fields(Dtest), 'M', Mlite, 'C', C, 'C_hom', C_hom, 'Cs33', Cs33);
    save(fullfile(dir2, 'results.mat'), '-struct', 's');

    res{ic} = struct('P2', P2, 'Dtest', Dtest, 'M', M, 'C', C);

end


%% ========================================================================
%  FIGURES (PNG), after the whole run
% ========================================================================

if P1.out.png
    fprintf('\n');
    fprintf('SAVING FIGURES ...\n');
    save_figures(P1, H, Dident, par1, res);
end

fprintf('\nDONE  (%.1f min)\n', toc(tAll)/60);

end


%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function H = run_homog(quick)
% homogenization of the periodic cell with all fields (stage_cache calls this only when the saved data are missing or outdated)
[C_hom, Cs33, info] = homogenize2D_PBC('case1', quick);
H = struct('C_hom', C_hom, 'Cs33', Cs33, 'info', info);
end


function fz = expand_force(force)
% stages to recompute: the forced ones and everything that depends on them
stages = {'homog', 'fullIdent', 'identification', 'fullTest', 'micro'};
down = struct('homog', {{'identification', 'micro'}}, ...
              'fullIdent', {{'identification', 'micro'}}, ...
              'identification', {{'micro'}}, ...
              'fullTest', {{}}, ...
              'micro', {{}});
fz = struct();
for i = 1:numel(stages), fz.(stages{i}) = false; end
for i = 1:numel(force)
    f = force{i};
    if strcmpi(f, 'all')
        for j = 1:numel(stages), fz.(stages{j}) = true; end
        continue
    end
    if ~isfield(fz, f)
        error('Unknown stage ''%s'' (use: %s, all).', f, strjoin(stages, ', '));
    end
    fz.(f) = true;
    dn = down.(f);
    for j = 1:numel(dn), fz.(dn{j}) = true; end
end
end


function D = strip_fields(D)
% remove the large mesh fields (they are saved in the fullscale folders)
if isfield(D.hom, 'mesh'), D.hom = rmfield(D.hom, 'mesh'); end
if isfield(D.het, 'mesh'), D.het = rmfield(D.het, {'mesh', 'isInc'}); end
end


function write_matrix(A, file)
folder = fileparts(file);
if ~isempty(folder) && exist(folder, 'dir') ~= 7, mkdir(folder); end
try
    writematrix(A, file);
catch
    dlmwrite(file, A, 'precision', '%.12g');
end
end


function write_comparison_csv(C, D, file)
% one row per (load, N): dominant edge displacement (full-scale, homogeneous, micromorphic) and all errors of compare_results
nN = numel(C.Ns);
nL = numel(C.names);
n = nN*nL;
loadName = cell(n, 1);
v = zeros(n, 14);
i = 0;
for l = 1:nL
    for iN = 1:nN
        i = i + 1;
        loadName{i} = C.names{l};
        v(i, :) = [C.Ns(iN), D.het(iN).Lcell, C.refDom(iN, l), C.homDom(iN, l), C.mmDom(iN, l), ...
                   C.errDisp(iN, l), C.errHomDisp(iN, l), C.errEn(iN, l), C.errHomEn(iN, l), ...
                   C.errPsi(iN, l), C.errSig(iN, l), C.errSigHom(iN, l), C.secF.micro(iN, l), C.secM.micro(iN, l)];
    end
end
names = {'N', 'cell', 'u_fullscale', 'u_homogeneous', 'u_micromorphic', 'err_disp_micro_pct', 'err_disp_hom_pct', 'err_energy_micro_pct', ...
         'err_energy_hom_pct', 'err_psi_pct', 'err_sigma_micro_pct', 'err_sigma_hom_pct', 'sec_force_micro_pct', 'sec_moment_micro_pct'};
T = [table(loadName, 'VariableNames', {'load'}), array2table(v, 'VariableNames', names)];
writetable(T, file);
end
