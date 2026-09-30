% ============================================================
%  main_identify_CbarChatG.m
%
%  Identification of the micromorphic stiffness matrices Cbar, Chat and the
%  curvature modulus G for the 2D matrix/inclusion RUC:
%
%    W = 1/2 psi'C psi + 1/2 gamma'Cbar gamma + psi'Chat gamma
%        + 1/2 (L/N)^2 kappa' diag(G11,G22,G22,G11,G33,G33) kappa
%
%    C = C_hom (known, PBC homogenization), psi = sym(psi), gamma = eps_M - psi,
%    kappa = grad psi, Voigt [a11 a22 2a12]
%    Cbar = [Cb11 Cb12 0; Cb12 Cb11 0; 0 0 Cb33]  positive definite
%    Chat = [Ch11 Ch12 0; Ch12 Ch11 0; 0 0 Ch33]  positive definite
%    G11 = G22 = G33 = G  (cfg.id.Gmode = 'single', 7 unknowns)
%    or G11, G22, G33     (cfg.id.Gmode = 'three',  9 unknowns)
%
%  from energy, integrated stress int sigma_ij and integrated double stress
%  int sigma_ij x_k of classical FE solves, separately for every N.
%
%    1) dataGeneration     FE solves for every load case and N:
%                          A1..A19, I1..I24, 2U_het, int sigma, int sigma x
%    2) smoothData         removes the numerical fluctuation over N
%                          (every quantity fitted on its own)
%    3) identifyCbarChatG  all N together: every parameter a monotone function
%                          p(N) = c0 + c1/N + c2/N^2 with c >= 0 (positive,
%                          non-increasing in N), Chat positive definite and the
%                          whole energy block [C Chat; Chat Cbar] positive
%                          definite at every N (convex problem, fmincon)
%
%  Files needed on the path: dataGeneration.m, smoothData.m,
%  identifyCbarChatG.m, plotAterms.m, homogenize2D_PBC.m (+ PDE Toolbox for meshing).
% ============================================================

clear; clc; close all;

%% ================= SETTINGS =================

cfg = struct();
cfg.outDir = fullfile(pwd, 'CbarChatG_out');
cfg.plots  = true;
cfg.plotAterms = true;     % after the fit: one PNG per load case with A1..A19 vs N
                           % (<outDir>/<region>/A_terms/<case>_A.png)

% ---- step 1: run the FE data generation, or reuse DataGeneration.mat ----
cfg.runDataGeneration = true;

% geometry (unit square specimen, N x N cells, circular inclusion per cell)
cfg.geom = struct('Lx_tot', 1.0, 'Ly_tot', 1.0, ...
                  'Rfrac', 6/19, ...        % inclusion radius / cell size
                  'meshFrac', 0.05, ...     % element size / cell size (RUC meshes)
                  'meshFracHom', 0.02);     % element size (Stage 1 and PBC homogenization)
cfg.Nlist = 1:8;                            % unit cells per side

% material, plane strain
cfg.mat = struct('E_m', 70000, 'nu_m', 0.33, 'E_i', 3500, 'nu_i', 0.33, 'thickness', 1.0);

% loading: raw 9-vectors [H11 H22 H12 G1_11 G1_22 G1_12 G2_11 G2_22 G2_12]
cfg.load  = struct('strain0', 0.01, 'gScaleFactor', 2);
cfg.cases = struct('source', 'library', ...       % 'library' | 'manual' | 'file'
                   'library', 'H3G16', ...        % H11, H22, H12 + 16 G cases
                   'manual', [0 0 0 0 0 0 1 0 0], ...
                   'numFromFile', 64, 'sampleFile', 'sample_points_d09.txt');

% translation averaging: mode 5 = nPeriod x nPeriod shifts over one whole cell
cfg.trans = struct('mode', 5, 'nPeriod', 8, 'marginFrac', 0.05, 'stepFrac', 0.05);

% fields
cfg.num.nGrid     = 121;       % evaluation grid points per side (minimum)
cfg.num.gridPerCell = 40;      % ... and at least this many per unit cell: nGrid(N) = max(121, 40N+1)
cfg.kappaMethod   = 'local';   % kappa = grad psi: 'local' = grad eps_M (exact, Stage-1 FE)
                               % - grad gamma (moving least squares, keeps boundary layers)
                               % | 'localPlain' (window on psi) | 'polyfit' | 'element' | 'gridFD'
cfg.kappaLocalWidth = 0.1;     % 'local': Gaussian width sigma = kappaLocalWidth * cell size;
                               % needs sigma >~ 0.6*cell/nPeriod to damp the shift ripple
                               % (0.1 with nPeriod = 8); larger values smooth macro gradients
cfg.kappaPolyDeg  = 2;         % 'polyfit': total degree

% regions: every integral also over the interior [a, L-a]^2, a = interiorBand * cell size
% (bulk without the KUBC boundary layer; empty for N <= 2*interiorBand)
cfg.interiorBand  = 0.5;
cfg.regions       = {'full', 'interior'};   % identified separately, results in subfolders

% ---- step 2: smoothing over N ----
cfg.smooth = struct('NminFit', 1, ...                    % N used from here on (1 = all)
                    'decayExponents', [0.5 1 1.5 2], ... % decaying quantities: sum c_k N^-p
                    'T4Exponents', [0 1 2], ...          % kappa squares: c0 + c1/N + c2/N^2, c >= 0
                    'freeExponents', [0 1 2]);           % all other quantities (free signs)

% ---- step 3: identification ----
cfg.id = struct('mode', 'monotone', ...          % 'monotone' (all N, p(N) decreasing, energy PD)% 'perN' | 'global' | 'both' (free in N)
                'trendExponents', [0 1 2], ...   % p(N) = sum c_k N^-e_k, c_k >= 0
                'energyPD', true, ...            % [C Chat; Chat Cbar] positive definite at every N
                'compareFree', true, ...         % also identify each N on its own, for comparison
                'Gmode', 'single', ...           % 'single': G11=G22=G33 (7 unknowns), 'three' (9)
                'wE', 1, 'wS', 1, 'wD', 1, ...   % channel weights (0 switches a channel off)
                'doubleStressModel', 'local+mu', ... % int sigma_FE x = int sigma_(gamma) x + int mu
                'tolInfo', 1e-2, ...             % skip groups that are ~0 (no information)
                'pdCbar', true, 'pdChat', true, ... % Cbar, Chat positive definite
                'requireEnergyPD', false, ...    % (modes perN/global) energy block PD as well
                'pdMargin', 1e-3, ...            % smallest eigenvalue >= pdMargin*max eig(C)
                'cases', []);                    % [] -> all load cases

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

% 2) + 3) smoothing and identification, for every region
Rall = struct();
for r = 1:numel(cfg.regions)
    reg  = cfg.regions{r};
    cfgR = cfg;  cfgR.region = reg;  cfgR.outDir = fullfile(cfg.outDir, reg);
    if ~isfolder(cfgR.outDir), mkdir(cfgR.outDir); end
    fprintf('\n######################## region: %s ########################\n', reg);
    try
        S = smoothData(D, cfgR);
        Rall.(reg) = identifyCbarChatG(S, cfgR);
        if cfg.plotAterms, plotAterms(S, cfgR); end
    catch err
        fprintf('[main] region ''%s'' skipped: %s\n', reg, err.message);
    end
end

%% ================= RESULT =================

fprintf('\n================ IDENTIFIED Cbar, Chat, G ================\n');
fprintf('  kappa: %s', cfg.kappaMethod);
if strcmpi(cfg.kappaMethod, 'local'), fprintf(' (sigma = %.2f cell)', cfg.kappaLocalWidth); end
fprintf(',  interior band = %.2f cell\n', cfg.interiorBand);
regs = fieldnames(Rall);
for r = 1:numel(regs)
    R = Rall.(regs{r});
    fprintf('\n  ---- region: %s ----\n', regs{r});
    fprintf('  C_hom = [%.6g %.6g 0; %.6g %.6g 0; 0 0 %.6g]\n', R.C(1,1), R.C(1,2), R.C(1,2), R.C(2,2), R.C(3,3));
    if isfield(R, 'perN')
        fprintf('  free per N (is each parameter N-independent?):\n');
        fprintf('  %4s', 'N');  fprintf(' %12s', R.paramNames{:});  fprintf('\n');
        for iN = 1:numel(R.Nlist)
            fprintf('  %4d', R.Nlist(iN));  fprintf(' %12.6g', R.perN.params(iN,:));  fprintf('\n');
        end
    end
    if isfield(R, 'monotone')
        fprintf('  monotone, energy positive definite:\n');
        fprintf('  %4s', 'N');  fprintf(' %12s', R.paramNames{:});  fprintf('\n');
        for iN = 1:numel(R.monotone.Nlist)
            fprintf('  %4d', R.monotone.Nlist(iN));  fprintf(' %12.6g', R.monotone.params(iN,:));  fprintf('\n');
        end
        fprintf('  p(N) = c0 + c1/N + c2/N^2 (rows c0, c1, c2):\n');
        for k = 1:numel(R.monotone.exponents)
            fprintf('  N^-%g', R.monotone.exponents(k));  fprintf(' %12.6g', R.monotone.coef(k,:));  fprintf('\n');
        end
        fprintf('  RMS misfit (energy / stress / double stress) at each N:\n');
        for iN = 1:numel(R.monotone.Nlist)
            fprintf('  %4d   %.4f  %.4f  %.4f\n', R.monotone.Nlist(iN), R.monotone.rmsPerChannel(iN,:));
        end
    end
    fprintf('  table: %s\n', fullfile(cfg.outDir, regs{r}, 'IdentifiedCbarChatG_vsN.csv'));
end
fprintf('==========================================================\n');
