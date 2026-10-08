function P = parameters(caseName, quick)
% =========================================================================
%  PARAMETERS - the ONE place with every input of the study (material, geometry, loads, test case, numerics, constraints, output).
%
%     P = parameters('case1')          P = parameters('case2')          P = parameters('case1', true)   (true = quick check)
%
%  Everything else in the code reads P; nothing is hard-coded elsewhere.  Edit the values below, or add a new case in section 9.
%  Units: m, N; stiffness in N/m^2 (plane strain, per unit thickness).  Strain / stress vectors are [11; 22; 12] with engineering shear.
% =========================================================================
if nargin < 1, caseName = 'case1'; end
if nargin < 2, quick = false; end
P.caseName = caseName;

%% 1. MATERIAL  (periodic composite: one circular inclusion in every square cell)
P.material.Em = 70000;   P.material.num = 0.33;      % matrix: Young's modulus, Poisson ratio
P.material.Ei = 3500;    P.material.nui = 0.33;      % inclusion: Young's modulus, Poisson ratio
P.material.radiusRatio = 6/19;                       % inclusion radius / cell edge
P.material.Chom = [];                                % Calculated by homogenize2D_PBC.m
P.material.thickness = 1.0;                          % plate thickness

%% 2. LOAD LIBRARY: traction on the edge x = Lx,  t(y) = [tx; ty] = sum_m coef(:,m+1) P_m(s),  s = 2y/Ly - 1  (Legendre P0, P1, P2)
P.loadLibrary.shear    = [0 0 0; -1 0 0];            % ty = -1               (uniform shear force)
P.loadLibrary.axial    = [1 0 0;  0 0 0];            % tx = +1               (uniform tension)
P.loadLibrary.bending  = [0 1 0;  0 0 0];            % tx = s                (end moment)
P.loadLibrary.shearpar = [0 0 0; -1 0 1];            % ty = -1.5 (1 - s^2)   (parabolic shear, unit mean)

%% 3. IDENTIFICATION PLATE  (plate clamped at x = 0, traction on x = Lx;  data of steps "identification")
P.ident.Lx = 1.0;   P.ident.Ly = 0.5;                % plate length and height
P.ident.NList = 1:6;                                 % N = number of cells across the height (cell size Ly/N); n_x = Lx*N/Ly must be an integer
P.ident.loads = {'shear', 'axial', 'bending'};       % load cases (names of the load library)
P.ident.psiMode = 'translation';                     % how psi (micro strain) is found from the heterogeneous plate:
                                                     %   'translation': mean strain over psiPeriod^2 translated positions of the inclusion lattice
                                                     %   'centred'    : strain of the plate with the inclusion centred in every cell (NO translation)
P.ident.psiPeriod = 12;                               % translations per direction (only for 'translation'): 4 -> 16 positions

%% 4. TEST PLATE  (another problem to test the identified parameters; ONLY the centred inclusion is solved, no translation)
P.test.Lx = 2.0;    P.test.Ly = 0.5;
P.test.NList = 1:4;
P.test.loads = {'shear', 'bending', 'shearpar'};
P.test.psiMode = 'centred';                          % keep 'centred': fields of the centred plate are compared with the micromorphic result
P.test.psiPeriod = 8;                                % (unused for 'centred')

%% 5. NUMERICS OF THE FINITE ELEMENT SOLUTIONS AND OF THE DATA
P.fe.meshFrac = 0.05;                                % element size of the heterogeneous plate = meshFrac * cell size
P.fe.NB = 10;                                        % bins per cell edge for the average of strain and stress
P.fe.psiDeg = [10 8];                                % Legendre degrees (x, y) of the smooth fields eps_macro, psi, sigma
P.fe.nyMesh = 30;                                    % micromorphic FE mesh of the test plate: elements along the height
P.fe.nyCalib = 15;                                   % micromorphic mesh (elements along the height) used in the refinement of the identification (many forward solves)
P.fe.nBinY = 30;                                     % bins along the height for the micromorphic fields and the comparison

%% 6. MICROMORPHIC MODEL AND CONSTRAINTS ON THE PARAMETERS
%  W = 1/2 psi'C_psi(N) psi + 1/2 gamma'C_gamma(N) gamma + psi'C_couple(N) gamma + 1/2 l^2 G |grad psi|^2,   gamma = eps_macro - psi
%  The three constitutive matrices have EXACTLY THE SAME parameterization:   C = [a b 0; b a 0; 0 0 d],   every entry  theta(N) = c0 + c1/N + c2/N^2
%  and can be chosen independently:   'free' = identified (9 coefficients),   'hom' = fixed to the homogenised stiffness C_hom.
P.model.Cpsi    = 'free';                            % C_psi(N):     'free' | 'hom'
P.model.Cgamma  = 'free';                            % C_gamma(N):   'free' | 'hom'
P.model.Ccouple = 'free';                            % C_couple(N):  'free' | 'hom'
P.model.ellFactor = 1.0;                             % internal length  l = ellFactor * cell size
P.model.boundFactor = 2;                             % every free C_ij(N) lies in  0 <= C_ij(N) <= boundFactor * C_hom,ij   (positive, below boundFactor*C_hom)
P.model.monotone = true;                             % every free C_ij(N) is monotone in N (increasing or decreasing, chosen by the fit)
P.model.limit = 'all';                               % 'all': C_ij(N -> infinity) = C_hom,ij for every free matrix  |  'none'
%  STABILITY AND CLOSENESS TO THE HOMOGENEOUS SOLUTION (linear constraints, mode by mode).  C_psi, C_gamma, C_couple = [a b 0; b a 0; 0 0 d] share the
%  modes  sum = a+b,  diff = a-b,  shear = d.  With A = C_psi, B = C_couple, C = C_gamma (scalars of one mode) the energy matrix is [A B; B C]:
%      0 <= B,   B <= A - margin*C_hom,   B <= C - margin*C_hom      =>  A C >= B^2: positive definite energy (stable forward problem)
%      long-wave stiffness  C_eff = C - (C-B)^2/(A+C-2B) = B + u v/(u+v)   (u = C-B, v = A-B)   satisfies   B <= C_eff <= min(A, C)
%  Hence for N >= Nhom:  C_couple >= (1-homTol)*C_hom  and  C_psi <= (1+homTol)*C_hom   keep C_eff within about homTol of C_hom (homogeneous limit).
P.model.stability = true;                            % impose the constraints above (false: may give an unstable energy)
P.model.margin = 0.01;                               % margin of the strict inequalities (in units of C_hom) for N in the identified range
P.model.homTol = 0.02;                               % closeness to the homogeneous solution for N >= Nhom
P.model.Nhom = 5;                                    % from this N on the micromorphic solution must be close to the homogeneous one
P.model.refine = true;                               % refine ALL parameters (coefficients and G) with the FORWARD micromorphic solution of the identification plate:
P.model.dispWeight = 1000;                           %   minimise  ||A c - b||^2/||b||^2 + dispWeight*mean(r_u^2 + r_E^2)  (r_u, r_E: relative errors of the edge
                                                     %   displacement and of the energy w.r.t. the full-scale plate), same constraints.  Larger = fits displacement better
P.model.refineIter = 40;                             % maximum Levenberg-Marquardt iterations
P.model.Grange = [-6 2];                             % G / C_hom,11 <= 10^Grange(2) in the refinement
P.model.G0 = 1000;                                   % reference value of G (only used to scale the equations)

%% 7. RESIDUALS (weak forms) USED IN THE IDENTIFICATION
P.residual.test = [6 4];                             % degrees (x, y) of the interior test functions
P.residual.nQuad = [40 28];                          % Gauss points (x, y) of the area integrals
P.residual.nEdge = 40;                               % Gauss points on the edges
P.residual.nTestEdge = 4;                            % highest Legendre polynomial used as test function on the edges
P.residual.sections = [0.25 0.5 0.75 1.0];           % sections x = s*Lx where force and moment equilibrium are imposed
P.residual.muEdges = [1 1 1 1];                      % edges where l^2 G grad(psi).n = 0 is imposed [x=Lx, x=0, y=Ly, y=0]
P.residual.weights = struct('macro', 1, 'micro', 1, 'trac', 1, 'mu', 1, 'energy', 100, 'force', 10, 'moment', 10);   % group weights

%% 8. OUTPUT
P.out.root = 'results';
if quick, P.out.root = 'results_quick'; end          % quick runs use their own folder: they never overwrite the saved data of a full run
P.out.folder = fullfile(P.out.root, caseName);       % everything of this case: fullscale/, identification/, micromorphic/, logs/, figures/, results.mat
P.out.homFolder = fullfile(P.out.root, 'homogenization');   % homogenization (C_hom, fields, figures): shared by all cases
P.out.reuse = true;                                  % reuse every saved stage whose input did not change (false = recompute everything)
P.out.version = 1;                                   % part of every cache key: increase it after changing the CODE (not the parameters) to discard all saved data
P.out.png = true;                                    % save PNG figures at the end of the run
P.out.pngDPI = 150;                                  % resolution of the PNG files
P.out.tolerance = 2;                                 % target of the displacement error (%)

%% 9. quick
if quick                                             % quick check of the installation
    P.ident.NList = 1:3;  P.ident.psiPeriod = 2;  P.test.NList = 1:3;  P.fe.meshFrac = 0.1;
end

%% 10. DERIVED (do not edit): traction coefficient arrays trc(:,:,l) = [tx; ty] of the load cases
P.ident.trc = load_matrix(P.loadLibrary, P.ident.loads);
P.test.trc  = load_matrix(P.loadLibrary, P.test.loads);
end

function trc = load_matrix(lib, names)
trc = zeros(2, 3, numel(names));
for l = 1:numel(names), trc(:, :, l) = lib.(names{l}); end
end
