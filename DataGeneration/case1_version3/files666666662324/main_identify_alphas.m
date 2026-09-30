% ============================================================
%  main_identify_alphas.m
%
%  Micromorphic parameter identification for the 2D matrix/inclusion RUC
%
%    W = am/2 psi'C psi + ah/2 gamma'C gamma + ac psi'C gamma
%        + al/2 (L/N)^2 kappa'A kappa
%
%  from energy, integrated stress and integrated double stress of classical
%  FE solves.
%
%    1) dataGeneration  FE solves for every load case and N (Stages 1-4):
%                       T1..T5, 2U_het, int sigma, int sigma_ij x_k (FE) and
%                       S1..S4, M1..M3 (integrals of C eps_M, C psi, C gamma, A kappa)
%    2) smoothData      removes the numerical fluctuation over N; every
%                       quantity fitted on its own (no relation imposed)
%    3) identifyAlphas  weighted least squares with a positive-definite
%                       energy constraint, separately for every N (from all
%                       load cases at that N) -> alpha_m(N), alpha_h(N),
%                       alpha_c(N), alpha_l(N), plus a trend c0 + c1/N + c2/N^2
%
%  Files needed on the path: dataGeneration.m, smoothData.m,
%  identifyAlphas.m, homogenize2D_PBC.m (+ PDE Toolbox for meshing).
% ============================================================

clear; clc; close all;

%% ================= SETTINGS =================

cfg = struct();
cfg.outDir = fullfile(pwd, 'Identification_out');
cfg.plots  = true;

% ---- step 1: run the FE data generation, or reuse DataGeneration.mat ----
cfg.runDataGeneration = true;

% geometry (unit square specimen, N x N cells, circular inclusion per cell)
cfg.geom = struct('Lx_tot', 1.0, 'Ly_tot', 1.0, ...
                  'Rfrac', 6/19, ...        % inclusion radius / cell size
                  'meshFrac', 0.05, ...     % element size / cell size (RUC meshes)
                  'meshFracHom', 0.02);     % element size (Stage 1 and PBC homogenization)
cfg.Nlist = 1:5;                            % unit cells per side

% material, plane strain
cfg.mat = struct('E_m', 70000, 'nu_m', 0.33, 'E_i', 3500, 'nu_i', 0.33, 'thickness', 1.0);

% loading: raw 9-vectors [H11 H22 H12 G1_11 G1_22 G1_12 G2_11 G2_22 G2_12]
cfg.load  = struct('strain0', 0.01, 'gScaleFactor', 2);
cfg.cases = struct('source', 'library', ...       % 'library' | 'manual' | 'file'
                   'library', 'H3BEND4', ...        % 'G6' | 'G16' | 'G21', prefix 'H3' adds H11, H22, H12
                   'manual', [0 0 0 0 0 0 1 0 0], ...
                   'numFromFile', 64, 'sampleFile', 'sample_points_d09.txt');

% translation averaging: mode 5 = nPeriod x nPeriod shifts over one whole cell
cfg.trans = struct('mode', 5, 'nPeriod', 8, 'marginFrac', 0.05, 'stepFrac', 0.05);

% fields
cfg.num.nGrid     = 121;       % evaluation grid points per side (odd)
cfg.kappaMethod   = 'polyfit'; % kappa = grad sym(psi): 'polyfit' | 'element' | 'gridFD'
cfg.kappaPolyDeg  = 2;
cfg.Akappa        = 'Chom';    % 'Chom': kappa'A kappa = kx'C kx + ky'C ky (alpha_l dimensionless)
                               % 'identity': A = I (alpha_l then has stiffness units)

% ---- step 2: smoothing over N ----
cfg.smooth = struct('NminFit', 2, ...                 % N used from here on
                    'decayExponents', [0.5 1 1.5 2], ... % monotone decay basis N^-p
                    'T4Exponents', [0 1 2], ...       % T4 = c0 + c1/N + c2/N^2
                    'freeExponents', [0 1 2], ...     % basis for the stress-integral vectors
                    'vectors', true);                 % smooth the stress integrals over N

% ---- step 3: identification ----
cfg.id = struct('mode', 'both', ...              % 'perN': alphas for every N  (alpha(N)) % 'global': one set for all N, 'both'
                'wE', 1, 'wS', 1, 'wD', 1, ...   % channel weights (0 switches a channel off)
                'stressModel', 'sigma_gamma', ... % int sigma_FE = int sigma_(gamma) = ah S3 + ac S2  % 'both': also = int s_(sym) = am S2 + ac S3
                'doubleStressModel', 'local+mu', ... % int sigma_FE x = ah M3 + ac M2 + al ell^2 S4  % 'mu': int sigma_FE x = al ell^2 S4 only
                'tolInfo', 1e-2, ...             % skip stress groups that are ~0 (no information)
                'fixAlphaM', true, 'alphaM', 1, ... % alpha_m = 1 (affine equilibrium ansatz)
                'pdMargin', 1e-3, ...             % am*ah - ac^2 >= pdMargin*am*ah
                'cases', []);                     % [] -> all load cases

%% ================= RUN =================

if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
dataFile = fullfile(cfg.outDir, 'DataGeneration.mat');

% 1) FE data
if cfg.runDataGeneration || ~isfile(dataFile)
    D = dataGeneration(cfg);
else
    D = load(dataFile);
    fprintf('[main] reusing %s\n', dataFile);
end

% 2) smoothing over N
S = smoothData(D, cfg);

% 3) identification
A = identifyAlphas(S, cfg);

%% ================= RESULT =================

fprintf('\n=============== IDENTIFIED PARAMETERS ===============\n');
fprintf('  (A_kappa = %s, stress = %s, double stress = %s, %d load cases)\n', ...
    cfg.Akappa, cfg.id.stressModel, cfg.id.doubleStressModel, numel(A.cases));
if isfield(A, 'perN')
    fprintf('  %4s %12s %12s %12s %12s\n', 'N', 'alpha_m', 'alpha_h', 'alpha_c', 'alpha_l');
    for iN = 1:numel(A.Nlist)
        fprintf('  %4d %12.6g %12.6g %12.6g %12.6g\n', A.Nlist(iN), A.perN.alphas(iN,:));
    end
    fprintf('  trend alpha(N) = %s:\n', A.trend.formula);
    for j = 1:4, fprintf('    %-8s = %s\n', A.names{j}, A.trend.text{j}); end
end
if isfield(A, 'global')
    fprintf('  one set for all N: am = %.6g, ah = %.6g, ac = %.6g, al = %.6g\n', A.global.alphas);
end
fprintf('  table: %s\n', fullfile(cfg.outDir, 'IdentifiedAlphas_vsN.csv'));
fprintf('=====================================================\n');
