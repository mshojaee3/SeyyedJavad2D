% ============================================================
%  RUC_FE_2D_LinearElastic_withHomRef.m
%  Standalone 2D LINEAR-ELASTIC FE solver: NX x NY RUC array, circular
%  inclusions, isotropic matrix + inclusion, KUBC via the SAME
%  u_i(X) = H_ij Xc_j + 1/2 G_ijk Xc_j Xc_k  polynomial used in
%  Main2D_LE.py / RUC_FE_2D.m.
%
%  NEW IN THIS VERSION: for each load case, in addition to the
%  heterogeneous (matrix+inclusion) KUBC solve, ONE extra FE solve is
%  run on a plain, SINGLE-MATERIAL domain using the homogenized
%  stiffness C_hom (and its out-of-plane row Cs33_hom) from
%  homogenize2D_PBC.m, under the SAME KUBC (H,G). This "homogeneous
%  reference" solution needs no translation ensemble (no microstructure
%  to translate) and, since the outer domain Lx_tot x Ly_tot is fixed
%  regardless of N and C_hom doesn't depend on N either, it is computed
%  ONCE (right after the load-case library is built) and reused
%  unchanged for every N in Nlist. Its results are written to a
%  SEPARATE file, Ensemble_Summary_NxN_homogenized.csv, next to each
%  N's Ensemble_Summary_NxN.csv (which is otherwise untouched) -- the
%  homogenized file has the SAME content for every N, since neither
%  translation nor RVE size affects a homogeneous solve.
%
%  LOAD CASES -- 9-VECTOR SAMPLED CONVENTION ONLY.
%  All named load cases (UniaxialX, PureShearSym, ... ) have been
%  removed. Every load case is now a 9-component raw vector
%       raw = [H11, H22, H12, G1_11, G1_22, G1_12, G2_11, G2_22, G2_12]
%  -- EXACTLY Main2D_LE.py's SAMPLED_PREFIX ('NC_k...') convention and
%  RUC_FE_2D.m's FALLBACK_EPS layout -- turned into physical H, G by:
%       [H11,H22,H12]                       = raw(1:3) * strain0
%       [G1_11,G1_22,G1_12,G2_11,G2_22,G2_12] = raw(4:9) * strain0 * gScaleFactor
%  H and G do NOT use the same scale factor: gScaleFactor (default 2)
%  was reverse-engineered by comparing against a reference Abaqus row
%  for this exact sample point -- with a naive uniform strain0 scaling,
%  Fbar/H matched Abaqus exactly but every recovered G component came
%  out exactly half of Abaqus's. See the gScaleFactor comment in the
%  USER INPUT section for the full diagnosis; re-verify it if your
%  reference data ever comes from a different driver.
%
%  Raw vectors come from EITHER of two sources (set loadCaseSource):
%    'manual' -- you supply one or more 9-vectors directly in
%                manualRawVectors below, e.g. a single pure-uniaxial
%                direction  [1,0,0,0,0,0,0,0,0]
%    'file'   -- the first numCasesFromFile rows of sample_points_d09.txt
%                (must sit in the SAME FOLDER as this .m file) are used.
%                That file is a hierarchically-ordered LOG-energy point
%                set on S^8, so its own header says "the first N rows
%                are themselves an optimized N-point set" -- i.e.
%                choosing numCasesFromFile = 64 automatically gives you
%                a good 64-point design, not an arbitrary truncation.
%  Every case is named 'NC_k%04d' (0-indexed), matching Main2D_LE.py's
%  SAMPLED_PREFIX naming so the two pipelines line up if ever compared.
%
%  TRANSLATION ENSEMBLE -- translationMode CONVENTION.
%  The inclusion-translation sweep used to reduce microstructure-
%  specific bias is now generated the way you specified, entirely in
%  units of the CELL size (Lx,Ly), so the sweep is the same shape for
%  any N:
%    translationMode = 1  diagonal sweep      zeta1(i) = zeta2(i)
%    translationMode = 2  full grid           zeta1(i) x zeta2(j)
%    translationMode = 3  "plus" cross        centre + steps along
%                                              +/-X and +/-Y separately
%    translationMode = 4  centred only        runList = [0,0], i.e. NO
%                                              translation sweep -- one
%                                              RUC, one solution
%  For modes 1-3 the final reported value per load case is the MEAN
%  (and, for ALLSE, also the STD) over every translation in runList;
%  for mode 4 it is trivially that single solution. (This ensemble
%  applies ONLY to the heterogeneous RUC solve -- the homogeneous
%  reference solve above has no microstructure and so no translation.)
%
%  ELEMENT: CPE6-equivalent 6-node quadratic triangle (T6), plane
%  strain, 3-point Gauss quadrature -- identical shape functions and
%  Gauss rule to t6ShapeDeriv()/gauss3() in RUC_FE_2D.m, and to
%  Main2D_LE.py's CPE6 Abaqus element.
%
%  CONSTITUTIVE LAW: small-strain isotropic linear elasticity, plane
%  strain, two phases (Matrix, Inclusion), each with its own (E,nu):
%       sigma = D(E,nu) * [e11; e22; gamma12],   gamma12 = 2*e12
%       sigma33 = nu*(sigma11+sigma22),          e33 = 0   (plane strain)
%  matching the Abaqus 'S'/'E' field-output convention used in
%  Main2D_LE.py's S_COMPS/E_COMPS = {*11,*22,*33,*12}.
%
%  OUTPUTS (one row per load case, averaged over the translation
%  ensemble -- same column meaning as Main2D_LE.py's per-job Summary):
%    ALLSE_mean/std          -- domain-integrated strain energy (RAW
%                               total, matching Abaqus's ALLSE history
%                               output -- NOT a density)
%    RF_LEFT_x..RF_TOP_y    -- per-side reaction-force sums, each side
%                               summed independently over its OWN node
%                               set (shared corners double-count, same
%                               convention as Main2D_LE.py)
%    Volume
%    Savg_11,Savg_22,Savg_33,Savg_12   -- Gauss-weighted volume averages
%    Eavg_11,Eavg_22,Eavg_33,Eavg_12
%    Uavg_1,Uavg_2           -- area(volume)-weighted mean-of-6-nodal-U
%                               per element (matches Main2D_LE.py's
%                               corrected, non-node-count Uavg)
%    ALLSE_hom, Savg_*_hom, Eavg_*_hom, Uavg_*_hom  -- NOT in this
%                               file: see Ensemble_Summary_NxN_homogenized.csv
%                               below, written alongside this file.
%  and, in a companion per-translation detail file:
%    Fbar_11..Fbar_22, G111..G222         -- Method A: boundary LSQ fit
%    Fbar_11_int..Fbar_22_int, G*_int     -- Method B: corrected
%                                            surface integral
%    Q111..Q222              -- double stress Q_ikl = <S_ik * Xc_l>_V
%                               via DOMAIN Gauss quadrature
%    fitResidMax, c1_fit, c2_fit          -- Method A diagnostics
%    H11..G2_12 (physical)   -- the actual H,G used for that row, for
%                               provenance/traceability
%
%  and, in a SEPARATE file per N, Ensemble_Summary_NxN_homogenized.csv
%  (same content for every N -- no translation/RVE-size dependence):
%  the SAME columns as the heterogeneous file above (Case, LoadCase,
%  strain0, StepTime, ALLSE, RF_*, Volume, nBndNodes, Savg/Eavg/Uavg,
%  Fbar/G by both methods, Q, fitResidMax/c1_fit/c2_fit, ALLSE_std
%  (trivially 0, single solve), and all diagnostics) EXCEPT NX, NY,
%  zeta1, zeta2, which don't apply to a homogeneous solve.
%
%  DEPENDENCY: MATLAB PDE Toolbox (createpde/geometryFromEdges/
%  generateMesh), used ONLY for mesh generation -- same as
%  RUC_FE_2D.m. If unavailable, replace buildRUCMesh() with your own
%  mesher returning `nodes` (Nx2) and `elems6` (Mx6, Abaqus/PDE-
%  Toolbox T6 node order: 3 corners, then 3 midsides on edges
%  1-2,2-3,3-1).
%
%  Also requires homogenize2D_PBC.m (same folder / MATLAB path) for
%  C_hom and Cs33_hom.
%
%  NOTE: written and reviewed carefully but NOT executed here (no
%  MATLAB/PDE Toolbox in this environment) -- run it and report back
%  any error; the most likely failure points are the decsg geometry
%  spec inside buildRUCMesh()/buildPlainMesh() and the
%  readmatrix('CommentStyle',...) call if your MATLAB version predates
%  it (a manual fallback parser is included for that case).
% ============================================================

clear; clc;

%% ============ USER INPUT ============

% ---- specimen size (FIXED regardless of NX,NY, same convention as
%      Main2D_LE.py / RUC_FE_2D.m) ----
Lx_tot = 1.0; Ly_tot = 1.0;

% ---- RVE sizes to sweep (each produces its own Ensemble_Summary file) ----
Nlist = [1,2,3];

% ---- geometry / mesh fractions (same names as Main2D_LE.py) ----
Rfrac    = 6/19;   % inclusion radius / unit-cell size
meshFrac = 0.05;   % element size / unit-cell size

% ---- loading amplitude (uniform scaling applied to every raw 9-vector,
%      see header comment above) ----
strain0 = 0.01;

% ---- G raw-to-physical scale factor -----------------------------------
% H is scaled as H = raw(1:3)*strain0. G is scaled as
% G = raw(4:9)*strain0*gScaleFactor. gScaleFactor=2 was reverse-engineered
% from a reference Abaqus run on this exact sample point: with
% gScaleFactor=1 (the naive "uniform strain0" scaling), Fbar (and hence
% H) matched the Abaqus row EXACTLY, but all 8 recovered G components
% were exactly half of Abaqus's -- a clean, exact factor of 2 across
% every component, not a rounding effect. This says the raw-to-physical
% mapping used by whatever produced that Abaqus reference (likely
% Run_LE_Sampled_Ensemble.m or an equivalent driver) scales the gradient
% part of the raw vector by 2*strain0, not strain0, even though the KUBC
% polynomial itself (u_i = H_ij Xc_j + 1/2 G_ijk Xc_j Xc_k) already
% carries its own explicit 1/2 on G -- i.e. the two are independent
% conventions and this script now matches the established one. If you
% ever regenerate the Abaqus reference or get it from a different
% driver, re-check this ratio (see the diagnosis in chat) before trusting
% gScaleFactor=2 blindly.
gScaleFactor = 2;

% ---- material: MATRIX and INCLUSION, plane strain, isotropic ----
E_m  = 70000.0;  nu_m = 0.33;
E_i  = 3500.0;  nu_i = 0.33;   % <-- set the real inclusion contrast here
THICKNESS = 1.0;

% ---- homogenized (effective) stiffness, via periodic BC on a single
%      unit cell (see homogenize2D_PBC.m) -- independent of N ----
mat.E_m = E_m; mat.nu_m = nu_m;
mat.E_i = E_i;  mat.nu_i = nu_i;
geom.Rfrac = Rfrac; geom.meshFrac = 0.02;

[C_hom, Cs33_hom] = homogenize2D_PBC(mat, geom);
fprintf('Homogenized C_hom (Voigt, plane strain, [S11;S22;S12]=C*[E11;E22;gamma12]) =\n');
disp(C_hom);
fprintf('Homogenized out-of-plane row Cs33_hom (Savg33 = Cs33_hom*[E11;E22;gamma12]) =\n');
disp(Cs33_hom);

% ------------------------------------------------------------
% LOAD CASE SETTINGS
% ------------------------------------------------------------
loadCaseSource = 'file';   % 'file' or 'manual'

% used when loadCaseSource = 'file': take the FIRST this many rows of
% sample_points_d09.txt (which must sit next to this .m file). Because
% that file is hierarchically ordered, this IS an optimized N-point set.
numCasesFromFile = 64;
sampleFileName   = 'sample_points_d09.txt';

% used when loadCaseSource = 'manual': one row per case, each a raw
% 9-vector [H11,H22,H12,G1_11,G1_22,G1_12,G2_11,G2_22,G2_12] on the
% UNSCALED (unit-sphere) scale -- multiplied by strain0 below, same as
% the 'file' source. Example: pure H11 stretch direction.
manualRawVectors = [ ...
    +5.72245728820545696e-01 -1.92916345731379196e-01 -1.25580798095818208e-01 -4.43331170730187241e-01 -2.71724749117376430e-01 -1.94577047958669619e-01 +5.34194488340966944e-01 -9.46217548967859950e-02 +1.30358930864310052e-01
];

% ------------------------------------------------------------
% TRANSLATION SETTINGS  (all relative to the CELL, so the grid is
% the same shape for any Nx, Ny)
% ------------------------------------------------------------
marginFrac      = 0.05;   % safety margin, fraction of cell size
stepFrac        = 0.05;   % translation step, fraction of cell size
translationMode =1;      % 1=diagonal, 2=full grid, 3=plus-cross, 4=centred only

% ------------------------------------------------------------
% PLOTTING  (manual mode only -- see header note below)
% ------------------------------------------------------------
% PLOT_RESULTS turns on, for every (N, translation, load case) actually
% solved WHILE loadCaseSource='manual', four figures:
%   (1) Sij field  -- S11,S22,S33,S12 contours on the mesh
%   (2) Eij field  -- E11,E22,E33,E12 contours on the mesh
%   (3) U field    -- U1,U2 contours on the mesh
%   (4) Reaction-force curves on LEFT/RIGHT/BOTTOM/TOP, i.e. the actual
%       nodal Rx(y),Ry(y) / Rx(x),Ry(x) distribution along each side --
%       not just their RF_*_x/y sums.
% This is meant for a single-shot manual run (one N, one translation --
% translationMode=4 -- one row in manualRawVectors); it is deliberately
% NOT enabled for loadCaseSource='file' since that is normally a
% 64-case x multi-translation sweep and would open hundreds of figures.
% If you do run a manual sweep over several N/translations/cases with
% plotting on, each figure is labelled with N, zeta, and the case name
% so they don't overwrite one another -- just expect one set per combo.
PLOT_RESULTS = false;
plotDeformScale = 0;   % 0 = plot on the undeformed mesh (recommended at
                        % strain0~1%); >0 = exaggerate displacement by
                        % this factor for a visibly deformed shape.

% ---- output ----
OUTDIR = fullfile(pwd, 'LE_Linear_Summary');
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

%% ============ BUILD THE LOAD-CASE LIBRARY (once, independent of N) ============

rawVecs = loadRawVectors(loadCaseSource, numCasesFromFile, manualRawVectors, sampleFileName);
nCases = size(rawVecs,1);
caseNames = arrayfun(@(k) sprintf('NC_k%04d', k-1), (1:nCases)', 'UniformOutput', false);

fprintf('Loaded %d load case(s) from source ''%s''.\n', nCases, loadCaseSource);

%% ============ HOMOGENEOUS (C_hom) REFERENCE SOLUTION ============
%  For each load case: ONE FE solve on a plain, single-material domain
%  (no inclusion) using C_hom/Cs33_hom, with the SAME KUBC (H,G) as the
%  heterogeneous RUC uses. Independent of the translation ensemble (no
%  microstructure to translate) AND of N (Lx_tot,Ly_tot are fixed
%  regardless of N, and C_hom doesn't depend on N either) -- so this
%  runs ONCE here and is reused, unchanged, for every N below. Every
%  column of the heterogeneous Ensemble_Summary table is reproduced
%  here too (RF_*, Fbar/G by both methods, Q, all diagnostics),
%  computed on the homogeneous mesh/solve -- ALLSE_std is trivially 0
%  (a single solve, no ensemble) and there is no NX/NY/zeta column
%  since none of those affect a homogeneous solve.

TOL_hom = 1e-6*max(Lx_tot,Ly_tot);
Xc_hom = 0.5*Lx_tot;  Yc_hom = 0.5*Ly_tot;
meshSizeHom = geom.meshFrac * min(Lx_tot,Ly_tot);

[nodesH, elems6H] = buildPlainMesh(Lx_tot, Ly_tot, meshSizeHom);
NnodesH = size(nodesH,1); NelemH = size(elems6H,1); ndofH = 2*NnodesH;

[GP_L1h, GP_L2h, GP_Wh] = gauss3();
elemDataH(NelemH) = struct('dNdX',[],'detJ0',[],'w',[],'Xgp',[],'Ygp',[]);
for e = 1:NelemH
    coords6 = nodesH(elems6H(e,:),:);
    elemDataH(e).dNdX = cell(1,3);
    elemDataH(e).detJ0 = zeros(1,3);
    elemDataH(e).w = GP_Wh;
    elemDataH(e).Xgp = zeros(1,3);
    elemDataH(e).Ygp = zeros(1,3);
    for g = 1:3
        [~,dNdX,detJ,Xg,Yg] = t6ShapeDeriv(GP_L1h(g),GP_L2h(g),coords6);
        if detJ <= 0
            error('Homogeneous mesh: non-positive Jacobian in element %d.', e);
        end
        elemDataH(e).dNdX{g} = dNdX;
        elemDataH(e).detJ0(g) = detJ;
        elemDataH(e).Xgp(g) = Xg;
        elemDataH(e).Ygp(g) = Yg;
    end
end

DelemH = repmat({C_hom}, NelemH, 1);
K_hom = assembleGlobalK_LE(elems6H, elemDataH, DelemH, THICKNESS, ndofH);
[leftNH,rightNH,botNH,topNH,bndAllH] = classifyBoundary(nodesH, Lx_tot, Ly_tot, TOL_hom);
V0_hom = Lx_tot*Ly_tot*THICKNESS;

for ic = 1:nCases
    [H_ic, G_ic] = assembleHG(rawVecs(ic,:), strain0, gScaleFactor);
    homRows(ic) = solveHomogeneousCase( ...
        K_hom, nodesH, elems6H, elemDataH, leftNH,rightNH,botNH,topNH,bndAllH, ...
        C_hom, Cs33_hom, THICKNESS, V0_hom, Lx_tot, Ly_tot, Xc_hom, Yc_hom, ndofH, ...
        H_ic, G_ic, caseNames{ic}, strain0, ic); %#ok<SAGROW>
end
fprintf('Homogeneous (C_hom) reference solved for %d load case(s) -- one FE solve each, reused for every N.\n', nCases);

% ---- homogenized reference table: same column names/order as the
%      heterogeneous Ensemble_Summary table (no _hom suffix -- the
%      file name marks it), minus NX/NY/zeta1/zeta2 (meaningless for a
%      homogeneous solve). Built once here; identical for every N.
Thom = struct2table(homRows);

%% ============ MAIN SWEEP: N -> translation ensemble -> load case ============

TOL_FRAC = 1e-6;

for N = Nlist

    Lx = Lx_tot / N;  Ly = Ly_tot / N;
    L_cell = min(Lx, Ly);
    R        = Rfrac    * L_cell;
    meshSize = meshFrac * L_cell;
    TOL = TOL_FRAC * max(Lx_tot, Ly_tot);
    Xc = 0.5*Lx_tot;  Yc = 0.5*Ly_tot;

    fprintf('=================================================================\n');
    fprintf(' RVE SIZE N = %d  (Lx=Ly=%.6f cells, L_cell=%.6f)\n', N, Lx, L_cell);
    fprintf('=================================================================\n');

    % ------------------------------------------------------------
    % translation settings  (all relative to the CELL, so the grid is
    % the same shape for any Nx, Ny)
    % ------------------------------------------------------------
    zeta1Critical = Lx/2 - R - marginFrac*Lx;
    zeta2Critical = Ly/2 - R - marginFrac*Ly;
    if zeta1Critical < 0 || zeta2Critical < 0
        error('Critical translation is negative: inclusion or margin too large for N=%d.', N);
    end
    stepX = stepFrac * Lx;
    stepY = stepFrac * Ly;
    zeta1Vals = -zeta1Critical : stepX : zeta1Critical;
    zeta2Vals = -zeta2Critical : stepY : zeta2Critical;

    runList = [];
    if translationMode == 1
        nDiag = min(numel(zeta1Vals), numel(zeta2Vals));
        for i = 1:nDiag
            runList = [runList; zeta1Vals(i), zeta2Vals(i)]; %#ok<AGROW>
        end
    elseif translationMode == 2
        for i = 1:numel(zeta1Vals)
            for j = 1:numel(zeta2Vals)
                runList = [runList; zeta1Vals(i), zeta2Vals(j)]; %#ok<AGROW>
            end
        end
    elseif translationMode == 3
        runList = [0, 0];
        for s = stepX:stepX:zeta1Critical
            runList = [runList; s, 0]; %#ok<AGROW>
        end
        for s = stepX:stepX:zeta1Critical
            runList = [runList; -s, 0]; %#ok<AGROW>
        end
        for s = stepY:stepY:zeta2Critical
            runList = [runList; 0, s]; %#ok<AGROW>
        end
        for s = stepY:stepY:zeta2Critical
            runList = [runList; 0, -s]; %#ok<AGROW>
        end
    elseif translationMode == 4
        % only centered inclusion -> one RUC = one solution
        runList = [0, 0];
    else
        error('Unknown translationMode.');
    end

    Ntrans = size(runList,1);
    fprintf('  translationMode = %d  ->  %d translation(s) per load case\n', translationMode, Ntrans);

    % ---- accumulate rows per load case across the translation sweep ----
    rowsByCase = cell(nCases,1);
    for ic = 1:nCases, rowsByCase{ic} = []; end

    for it = 1:Ntrans

        zeta1 = runList(it,1);
        zeta2 = runList(it,2);

        fprintf('  -- translation %d/%d: zeta1=%.6f zeta2=%.6f --\n', it, Ntrans, zeta1, zeta2);

        % ---- build ONE mesh / stiffness matrix for this (N, zeta) pair;
        %      reused for every load case (only the BC changes) ----
        centers = buildInclusionCenters(N, N, Lx, Ly, zeta1, zeta2, Lx_tot, Ly_tot, R);
        [nodes, elems6, isInclusionElem] = buildRUCMesh(Lx_tot, Ly_tot, centers, R, meshSize);
        Nnodes = size(nodes,1); Nelem = size(elems6,1);

        [GP_L1, GP_L2, GP_W] = gauss3();
        clear elemData;
        elemData(Nelem) = struct('dNdX',[],'detJ0',[],'w',[],'Xgp',[],'Ygp',[]);
        for e = 1:Nelem
            coords6 = nodes(elems6(e,:),:);
            elemData(e).dNdX = cell(1,3);
            elemData(e).detJ0 = zeros(1,3);
            elemData(e).w = GP_W;
            elemData(e).Xgp = zeros(1,3);
            elemData(e).Ygp = zeros(1,3);
            for g = 1:3
                [~,dNdX,detJ,Xg,Yg] = t6ShapeDeriv(GP_L1(g),GP_L2(g),coords6);
                if detJ <= 0
                    error('Element %d has non-positive reference Jacobian at GP %d (detJ=%.4g).', e, g, detJ);
                end
                elemData(e).dNdX{g} = dNdX;
                elemData(e).detJ0(g) = detJ;
                elemData(e).Xgp(g) = Xg;
                elemData(e).Ygp(g) = Yg;
            end
        end

        [leftN, rightN, botN, topN, bndAll] = classifyBoundary(nodes, Lx_tot, Ly_tot, TOL);

        D_m = planeStrainD(E_m, nu_m);
        D_i = planeStrainD(E_i, nu_i);
        nu_e  = zeros(Nelem,1);
        Delem = cell(Nelem,1);
        for e = 1:Nelem
            if isInclusionElem(e)
                Delem{e} = D_i; nu_e(e) = nu_i;
            else
                Delem{e} = D_m; nu_e(e) = nu_m;
            end
        end

        ndof = 2*Nnodes;
        K = assembleGlobalK_LE(elems6, elemData, Delem, THICKNESS, ndof);
        V0 = Lx_tot * Ly_tot * THICKNESS;

        for ic = 1:nCases

            [H, G] = assembleHG(rawVecs(ic,:), strain0, gScaleFactor);

            [isDirichlet, uBC] = buildKUBC_BC(H, G, nodes, bndAll, Xc, Yc, ndof);

            % ---- linear solve (Dirichlet elimination) ----
            freeDofs = find(~isDirichlet);
            dirDofs  = find(isDirichlet);
            u = zeros(ndof,1);
            u(dirDofs) = uBC(dirDofs);
            if ~isempty(freeDofs)
                Kff = K(freeDofs, freeDofs);
                Kfd = K(freeDofs, dirDofs);
                rhs = -Kfd * u(dirDofs);
                u(freeDofs) = Kff \ rhs;
            end
            Rfull = K * u;

            [RF_L_x, RF_L_y] = sideReactionForces(Rfull, leftN);
            [RF_R_x, RF_R_y] = sideReactionForces(Rfull, rightN);
            [RF_B_x, RF_B_y] = sideReactionForces(Rfull, botN);
            [RF_T_x, RF_T_y] = sideReactionForces(Rfull, topN);

            [Savg, Eavg, ALLSE, Q1_S, Q2_S] = domainAverages_LE( ...
                u, elems6, elemData, Delem, nu_e, THICKNESS, V0, Xc, Yc);

            Q111 = Q1_S(1); Q112 = Q2_S(1);   % from S11
            Q221 = Q1_S(2); Q222 = Q2_S(2);   % from S22
            Q121 = Q1_S(3); Q122 = Q2_S(3);   % from S12
            Q211 = Q1_S(3); Q212 = Q2_S(3);   % from S21 = S12 (symmetric)

            Uavg = elementAveragedU(u, elems6, elemData, THICKNESS);

            [Fbar, Gm, fitResidMax, c1_fit, c2_fit] = boundaryFbarG_LSQ( ...
                u, nodes, bndAll, Xc, Yc);
            [Fbar_int, G_int] = boundaryFbarG_Integral( ...
                u, nodes, leftN, rightN, botN, topN, Lx_tot, Ly_tot, Xc, Yc);

            % ---- diagnostics, same Gauss-quadrature / boundary-reaction
            %      cross-check philosophy as RUC_FE_2D.m's postprocessFrame,
            %      ported to the linear-elastic (Cauchy stress) setting ----
            [HillResid, QSymResid, ForceBalanceResid, KubcErrF, KubcErrG, EnergyConsistencyResid] = ...
                computeDiagnostics_LE(Rfull, nodes, bndAll, Xc, Yc, Savg, Q1_S, Q2_S, V0, ...
                    Fbar, Gm, H, G, u, K, ALLSE);

            row = struct();
            row.LoadCase = caseNames{ic};
            row.zeta1 = zeta1; row.zeta2 = zeta2;
            row.ALLSE = ALLSE;
            row.RF_LEFT_x = RF_L_x;   row.RF_LEFT_y = RF_L_y;
            row.RF_RIGHT_x = RF_R_x;  row.RF_RIGHT_y = RF_R_y;
            row.RF_BOTTOM_x = RF_B_x; row.RF_BOTTOM_y = RF_B_y;
            row.RF_TOP_x = RF_T_x;    row.RF_TOP_y = RF_T_y;
            row.Volume = V0;
            row.nBndNodes = numel(bndAll);
            row.Savg_11 = Savg(1); row.Savg_22 = Savg(2);
            row.Savg_33 = Savg(3); row.Savg_12 = Savg(4);
            row.Eavg_11 = Eavg(1); row.Eavg_22 = Eavg(2);
            row.Eavg_33 = Eavg(3); row.Eavg_12 = Eavg(4);
            row.Uavg_1 = Uavg(1);  row.Uavg_2 = Uavg(2);
            row.Fbar_11 = Fbar(1,1); row.Fbar_12 = Fbar(1,2);
            row.Fbar_21 = Fbar(2,1); row.Fbar_22 = Fbar(2,2);
            row.G111 = Gm(1,1,1); row.G112 = Gm(1,1,2);
            row.G121 = Gm(1,2,1); row.G122 = Gm(1,2,2);
            row.G211 = Gm(2,1,1); row.G212 = Gm(2,1,2);
            row.G221 = Gm(2,2,1); row.G222 = Gm(2,2,2);
            row.Fbar_11_int = Fbar_int(1,1); row.Fbar_12_int = Fbar_int(1,2);
            row.Fbar_21_int = Fbar_int(2,1); row.Fbar_22_int = Fbar_int(2,2);
            row.G111_int = G_int(1,1,1); row.G112_int = G_int(1,1,2);
            row.G121_int = G_int(1,2,1); row.G122_int = G_int(1,2,2);
            row.G211_int = G_int(2,1,1); row.G212_int = G_int(2,1,2);
            row.G221_int = G_int(2,2,1); row.G222_int = G_int(2,2,2);
            row.Q111 = Q111; row.Q112 = Q112; row.Q121 = Q121; row.Q122 = Q122;
            row.Q211 = Q211; row.Q212 = Q212; row.Q221 = Q221; row.Q222 = Q222;
            row.fitResidMax = fitResidMax; row.c1_fit = c1_fit; row.c2_fit = c2_fit;
            row.H11 = H(1,1); row.H22 = H(2,2); row.H12 = H(1,2);
            row.G1_11 = G(1,1,1); row.G1_22 = G(1,2,2); row.G1_12 = G(1,1,2);
            row.G2_11 = G(2,1,1); row.G2_22 = G(2,2,2); row.G2_12 = G(2,1,2);
            row.HillResid = HillResid; row.QSymResid = QSymResid;
            row.ForceBalanceResid = ForceBalanceResid;
            row.KubcErrF = KubcErrF; row.KubcErrG = KubcErrG;
            row.EnergyConsistencyResid = EnergyConsistencyResid;

            rowsByCase{ic} = [rowsByCase{ic}; row]; %#ok<AGROW>

            % ---- optional field / reaction-force plots (manual mode only) ----
            if PLOT_RESULTS && strcmpi(loadCaseSource, 'manual')
                tagStr = sprintf('N=%d  case=%s  zeta=(%.4f,%.4f)', ...
                    N, caseNames{ic}, zeta1, zeta2);
                [Sfield, Efield, Ufield] = computeNodalFields(u, nodes, elems6, Delem, nu_e);
                plotStressStrainDisp(nodes, elems6, Sfield, Efield, Ufield, ...
                    plotDeformScale, tagStr);
                plotReactionForceCurves(Rfull, nodes, leftN, rightN, botN, topN, tagStr);
            end
        end
    end

    % ---- ENSEMBLE AVERAGE over translations, per load case ----
    % Every physical quantity is the MEAN over the translation ensemble;
    % every residual/diagnostic quantity is the MAX over the ensemble
    % (a residual is meant to flag a problem -- averaging it away would
    % hide an isolated bad translation). ALLSE_std is reported alongside
    % as its own diagnostic: a large spread across translations signals
    % the RVE is too small / under-sampled for that load case.
    outCase        = (1:nCases)';
    outLoadCase    = caseNames;
    outNX          = N*ones(nCases,1);
    outNY          = N*ones(nCases,1);
    outStrain0     = strain0*ones(nCases,1);
    outStepTime    = ones(nCases,1);   % single linear step -> frame time 1.0
    outzeta1       = zeros(nCases,1);
    outzeta2       = zeros(nCases,1);
    outALLSE       = zeros(nCases,1);
    outALLSE_std   = zeros(nCases,1);
    outRF          = zeros(nCases,8);
    outVolume      = zeros(nCases,1);
    outnBnd        = zeros(nCases,1);
    outSavg        = zeros(nCases,4);
    outEavg        = zeros(nCases,4);
    outUavg        = zeros(nCases,2);
    outFbar        = zeros(nCases,4);
    outG           = zeros(nCases,8);
    outFbar_int    = zeros(nCases,4);
    outG_int       = zeros(nCases,8);
    outQ           = zeros(nCases,8);
    outFitResidMax = zeros(nCases,1);
    outC1fit       = zeros(nCases,1);
    outC2fit       = zeros(nCases,1);
    outHillResid   = zeros(nCases,1);
    outQSymResid   = zeros(nCases,1);
    outForceBalanceResid = zeros(nCases,1);
    outKubcErrF    = zeros(nCases,1);
    outKubcErrG    = zeros(nCases,1);
    outEnergyResid = zeros(nCases,1);

    for ic = 1:nCases
        R_ = rowsByCase{ic};

        outzeta1(ic) = mean([R_.zeta1]);
        outzeta2(ic) = mean([R_.zeta2]);

        allse = [R_.ALLSE];
        outALLSE(ic)     = mean(allse);
        outALLSE_std(ic) = std(allse);

        outRF(ic,:) = mean([ [R_.RF_LEFT_x]', [R_.RF_LEFT_y]', ...
                              [R_.RF_RIGHT_x]',[R_.RF_RIGHT_y]', ...
                              [R_.RF_BOTTOM_x]',[R_.RF_BOTTOM_y]', ...
                              [R_.RF_TOP_x]',   [R_.RF_TOP_y]' ], 1);

        outVolume(ic) = mean([R_.Volume]);
        outnBnd(ic)   = round(mean([R_.nBndNodes]));

        outSavg(ic,:) = mean([ [R_.Savg_11]',[R_.Savg_22]',[R_.Savg_33]',[R_.Savg_12]' ], 1);
        outEavg(ic,:) = mean([ [R_.Eavg_11]',[R_.Eavg_22]',[R_.Eavg_33]',[R_.Eavg_12]' ], 1);
        outUavg(ic,:) = mean([ [R_.Uavg_1]', [R_.Uavg_2]' ], 1);

        outFbar(ic,:) = mean([ [R_.Fbar_11]',[R_.Fbar_12]',[R_.Fbar_21]',[R_.Fbar_22]' ], 1);
        outG(ic,:)    = mean([ [R_.G111]',[R_.G112]',[R_.G121]',[R_.G122]', ...
                                [R_.G211]',[R_.G212]',[R_.G221]',[R_.G222]' ], 1);

        outFbar_int(ic,:) = mean([ [R_.Fbar_11_int]',[R_.Fbar_12_int]', ...
                                    [R_.Fbar_21_int]',[R_.Fbar_22_int]' ], 1);
        outG_int(ic,:)    = mean([ [R_.G111_int]',[R_.G112_int]',[R_.G121_int]',[R_.G122_int]', ...
                                    [R_.G211_int]',[R_.G212_int]',[R_.G221_int]',[R_.G222_int]' ], 1);

        outQ(ic,:) = mean([ [R_.Q111]',[R_.Q112]',[R_.Q121]',[R_.Q122]', ...
                             [R_.Q211]',[R_.Q212]',[R_.Q221]',[R_.Q222]' ], 1);

        outFitResidMax(ic) = max([R_.fitResidMax]);
        outC1fit(ic) = mean([R_.c1_fit]);
        outC2fit(ic) = mean([R_.c2_fit]);

        outHillResid(ic)         = max([R_.HillResid]);
        outQSymResid(ic)         = max([R_.QSymResid]);
        outForceBalanceResid(ic) = max([R_.ForceBalanceResid]);
        outKubcErrF(ic)          = max([R_.KubcErrF]);
        outKubcErrG(ic)          = max([R_.KubcErrG]);
        outEnergyResid(ic)       = max([R_.EnergyConsistencyResid]);

        % ---- final per-case results, printed to the MATLAB console ----
        fprintf('---------------------------------------------------------------\n');
        fprintf(' Case %d  [%s]   N=%dx%d   (mean over %d translation(s))\n', ...
            ic, caseNames{ic}, N, N, Ntrans);
        fprintf('---------------------------------------------------------------\n');
        fprintf('  ALLSE           = %.6e  (std over translations %.2e)   ALLSE_hom = %.6e\n', ...
            outALLSE(ic), outALLSE_std(ic), Thom.ALLSE(ic));
        fprintf('  Savg [11,22,33,12]     = [%.4e  %.4e  %.4e  %.4e]\n', outSavg(ic,:));
        fprintf('  Savg_hom [11,22,33,12] = [%.4e  %.4e  %.4e  %.4e]\n', ...
            Thom.Savg_11(ic), Thom.Savg_22(ic), Thom.Savg_33(ic), Thom.Savg_12(ic));
        fprintf('  Eavg [11,22,33,12] = [%.4e  %.4e  %.4e  %.4e]\n', outEavg(ic,:));
        fprintf('  Uavg [1,2]         = [%.4e  %.4e]\n', outUavg(ic,:));
        fprintf('  RF   L=(%.3e,%.3e)  R=(%.3e,%.3e)  B=(%.3e,%.3e)  T=(%.3e,%.3e)\n', outRF(ic,:));
        fprintf('  Fbar (LSQ) = [%.6f %.6f; %.6f %.6f]   Fbar (int) = [%.6f %.6f; %.6f %.6f]\n', ...
            outFbar(ic,:), outFbar_int(ic,:));
        fprintf('  G (LSQ)    = [%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e]\n', outG(ic,:));
        fprintf('  Q          = [%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e]\n', outQ(ic,:));
        fprintf(['  Diagnostics: HillResid=%.2e  QSymResid=%.2e  ForceBalanceResid=%.2e\n' ...
                 '               KubcErrF=%.2e   KubcErrG=%.2e   EnergyConsistencyResid=%.2e' ...
                 '   fitResidMax=%.2e\n'], ...
            outHillResid(ic), outQSymResid(ic), outForceBalanceResid(ic), ...
            outKubcErrF(ic), outKubcErrG(ic), outEnergyResid(ic), outFitResidMax(ic));
    end

    T = table( outCase, outLoadCase, outNX, outNY, outzeta1, outzeta2, outStrain0, outStepTime, ...
        outALLSE, ...
        outRF(:,1),outRF(:,2),outRF(:,3),outRF(:,4), ...
        outRF(:,5),outRF(:,6),outRF(:,7),outRF(:,8), ...
        outVolume, outnBnd, ...
        outSavg(:,1),outSavg(:,2),outSavg(:,3),outSavg(:,4), ...
        outEavg(:,1),outEavg(:,2),outEavg(:,3),outEavg(:,4), ...
        outUavg(:,1),outUavg(:,2), ...
        outFbar(:,1),outFbar(:,2),outFbar(:,3),outFbar(:,4), ...
        outG(:,1),outG(:,2),outG(:,3),outG(:,4),outG(:,5),outG(:,6),outG(:,7),outG(:,8), ...
        outFbar_int(:,1),outFbar_int(:,2),outFbar_int(:,3),outFbar_int(:,4), ...
        outG_int(:,1),outG_int(:,2),outG_int(:,3),outG_int(:,4), ...
        outG_int(:,5),outG_int(:,6),outG_int(:,7),outG_int(:,8), ...
        outQ(:,1),outQ(:,2),outQ(:,3),outQ(:,4),outQ(:,5),outQ(:,6),outQ(:,7),outQ(:,8), ...
        outFitResidMax, outC1fit, outC2fit, ...
        outALLSE_std, outHillResid, outQSymResid, outForceBalanceResid, ...
        outKubcErrF, outKubcErrG, outEnergyResid, ...
        'VariableNames', { ...
        'Case','LoadCase','NX','NY','zeta1','zeta2','strain0','StepTime', ...
        'ALLSE', ...
        'RF_LEFT_x','RF_LEFT_y','RF_RIGHT_x','RF_RIGHT_y', ...
        'RF_BOTTOM_x','RF_BOTTOM_y','RF_TOP_x','RF_TOP_y', ...
        'Volume','nBndNodes', ...
        'Savg_11','Savg_22','Savg_33','Savg_12', ...
        'Eavg_11','Eavg_22','Eavg_33','Eavg_12', ...
        'Uavg_1','Uavg_2', ...
        'Fbar_11','Fbar_12','Fbar_21','Fbar_22', ...
        'G111','G112','G121','G122','G211','G212','G221','G222', ...
        'Fbar_11_int','Fbar_12_int','Fbar_21_int','Fbar_22_int', ...
        'G111_int','G112_int','G121_int','G122_int', ...
        'G211_int','G212_int','G221_int','G222_int', ...
        'Q111','Q112','Q121','Q122','Q211','Q212','Q221','Q222', ...
        'fitResidMax','c1_fit','c2_fit', ...
        'ALLSE_std','HillResid','QSymResid','ForceBalanceResid', ...
        'KubcErrF','KubcErrG','EnergyConsistencyResid'});

    outName = fullfile(OUTDIR, sprintf('Ensemble_Summary_%dX%d.csv', N, N));
    writetable(T, outName);
    fprintf('Wrote %s\n', outName);

    % ---- separate homogenized-reference file for this N -- SAME
    %      content as for every other N (Thom was built once above,
    %      before the N loop); written again here only so it sits
    %      alongside its matching Ensemble_Summary_NxN.csv.
    homOutName = fullfile(OUTDIR, sprintf('Ensemble_Summary_%dX%d_homogenized.csv', N, N));
    writetable(Thom, homOutName);
    fprintf('Wrote %s\n', homOutName);

    % ---- full per-translation detail (Fbar/G/Q/fit diagnostics + the
    %      physical H,G actually used) for inspection ----
    detailRows = vertcat(rowsByCase{:});
    Tdetail = struct2table(detailRows);
    detailName = fullfile(OUTDIR, sprintf('LoadCaseDetail_%dX%d.csv', N, N));
    writetable(Tdetail, detailName);
    fprintf('Wrote %s\n', detailName);

end

fprintf('\nAll RVE sizes completed.\n');


% ============================================================
% LOAD-CASE VECTOR LOADING  (manual list, or head of the sample file)
% ============================================================

function rawVecs = loadRawVectors(source, nCases, manualVecs, fileName)
    if strcmpi(source, 'manual')
        rawVecs = manualVecs;
        if size(rawVecs,2) ~= 9
            error(['manualRawVectors must have 9 columns ' ...
                   '[H11,H22,H12,G1_11,G1_22,G1_12,G2_11,G2_22,G2_12].']);
        end
    elseif strcmpi(source, 'file')
        scriptDir = fileparts(mfilename('fullpath'));
        filePath = fullfile(scriptDir, fileName);
        if ~isfile(filePath)
            error('Sample points file not found next to this script: %s', filePath);
        end
        allVecs = readSamplePoints(filePath);
        if nCases > size(allVecs,1)
            error('Requested %d cases but %s only has %d data rows.', ...
                nCases, fileName, size(allVecs,1));
        end
        % the file is hierarchically ordered so its first N rows are
        % themselves an optimized N-point set -- simply take the head.
        rawVecs = allVecs(1:nCases, :);
    else
        error('Unknown loadCaseSource: %s (use ''file'' or ''manual'').', source);
    end
end

function vecs = readSamplePoints(filePath)
    try
        vecs = readmatrix(filePath, 'CommentStyle', '#');
    catch
        % fallback for MATLAB versions where readmatrix lacks
        % 'CommentStyle' -- parse by hand, skipping '#' header lines.
        fid = fopen(filePath, 'rt');
        if fid < 0
            error('Could not open %s.', filePath);
        end
        rows = {};
        tline = fgetl(fid);
        while ischar(tline)
            s = strtrim(tline);
            if ~isempty(s) && s(1) ~= '#'
                rows{end+1} = sscanf(s, '%f')'; %#ok<AGROW>
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        vecs = cell2mat(rows(:));
    end
    if size(vecs,2) ~= 9
        error('Expected 9 columns in %s, found %d.', filePath, size(vecs,2));
    end
end


% ============================================================
% RAW 9-VECTOR -> PHYSICAL (H, G)   (same layout as Main2D_LE.py's
% SAMPLED_PREFIX branch / RUC_FE_2D.m's FALLBACK_EPS). H uses strain0;
% G uses strain0*gScaleFactor -- see the gScaleFactor comment in the
% USER INPUT section for why these are NOT the same scale.
% ============================================================

function [H, G] = assembleHG(raw9, strain0, gScaleFactor)
    vH = raw9(1:3) * strain0;
    vG = raw9(4:9) * (strain0 * gScaleFactor);
    H = zeros(2,2);
    H(1,1) = vH(1); H(2,2) = vH(2); H(1,2) = vH(3); H(2,1) = vH(3);
    G = zeros(2,2,2);   % G(i,j,k), symmetric in (j,k)
    G(1,1,1) = vG(1); G(1,2,2) = vG(2); G(1,1,2) = vG(3); G(1,2,1) = vG(3);
    G(2,1,1) = vG(4); G(2,2,2) = vG(5); G(2,1,2) = vG(6); G(2,2,1) = vG(6);
end


% ============================================================
% BOUNDARY CONDITION ASSEMBLY  (full KUBC on every boundary node --
% the only case left, now that all named/mixed cases are removed)
% ============================================================

function u = polyDisp(H, G, X, Y, Xc, Yc)
    xr = X - Xc; yr = Y - Yc;
    u = zeros(2,1);
    for i = 1:2
        u(i) = H(i,1)*xr + H(i,2)*yr ...
             + 0.5*G(i,1,1)*xr^2 + G(i,1,2)*xr*yr + 0.5*G(i,2,2)*yr^2;
    end
end

function [isDirichlet, uBC] = buildKUBC_BC(H, G, nodes, bndAll, Xc, Yc, ndof)
    isDirichlet = false(ndof,1);
    uBC = zeros(ndof,1);
    for k = 1:numel(bndAll)
        n = bndAll(k);
        uval = polyDisp(H, G, nodes(n,1), nodes(n,2), Xc, Yc);
        isDirichlet(2*n-1) = true; uBC(2*n-1) = uval(1);
        isDirichlet(2*n)   = true; uBC(2*n)   = uval(2);
    end
end


% ============================================================
% LINEAR ELASTIC CONSTITUTIVE MATRIX (plane strain)
% ============================================================

function D = planeStrainD(E, nu)
    % strain order [e11; e22; gamma12], gamma12 = 2*e12
    c = E / ((1+nu)*(1-2*nu));
    D = c * [ 1-nu,  nu,        0; ...
              nu,    1-nu,      0; ...
              0,     0,   (1-2*nu)/2 ];
end


% ============================================================
% ELEMENT / GLOBAL STIFFNESS  (linear, assembled ONCE per mesh)
% ============================================================

function Ke = elementStiffnessLE(ed, Dmat, THICKNESS)
    Ke = zeros(12,12);
    for g = 1:3
        dNdX = ed.dNdX{g};
        B = zeros(3,12);
        for a = 1:6
            B(1, 2*a-1) = dNdX(a,1);
            B(2, 2*a)   = dNdX(a,2);
            B(3, 2*a-1) = dNdX(a,2);
            B(3, 2*a)   = dNdX(a,1);
        end
        wdv = ed.w(g) * ed.detJ0(g) * THICKNESS;
        Ke = Ke + wdv * (B.' * Dmat * B);
    end
end

function K = assembleGlobalK_LE(elems6, elemData, Delem, THICKNESS, ndof)
    Nelem = size(elems6,1);
    Ii = zeros(Nelem*144,1); Jj = zeros(Nelem*144,1); Vv = zeros(Nelem*144,1);
    ptr = 0;
    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        Ke = elementStiffnessLE(elemData(e), Delem{e}, THICKNESS);
        for a = 1:12
            for b = 1:12
                ptr = ptr+1;
                Ii(ptr) = gdofs(a); Jj(ptr) = gdofs(b); Vv(ptr) = Ke(a,b);
            end
        end
    end
    K = sparse(Ii(1:ptr), Jj(1:ptr), Vv(1:ptr), ndof, ndof);
end


% ============================================================
% DOMAIN AVERAGES: Savg, Eavg, ALLSE (total), Q via Gauss quadrature
% ============================================================

function [Savg, Eavg, ALLSE, Q1_S, Q2_S] = domainAverages_LE( ...
        u, elems6, elemData, Delem, nu_e, THICKNESS, V0, Xc, Yc)

    Nelem = size(elems6,1);
    Sacc = zeros(4,1);   % S11,S22,S33,S12
    Eacc = zeros(4,1);   % E11,E22,E33(=0),E12
    ALLSE = 0.0;

    Q1acc = zeros(3,1);  % moments against Xc, for S11,S22,S12
    Q2acc = zeros(3,1);  % moments against Yc, for S11,S22,S12

    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        Dmat = Delem{e};
        nu = nu_e(e);

        for g = 1:3
            dNdX = elemData(e).dNdX{g};
            B = zeros(3,12);
            for a = 1:6
                B(1, 2*a-1) = dNdX(a,1);
                B(2, 2*a)   = dNdX(a,2);
                B(3, 2*a-1) = dNdX(a,2);
                B(3, 2*a)   = dNdX(a,1);
            end
            eeng = B * ue;                 % [e11; e22; gamma12]
            s = Dmat * eeng;               % [s11; s22; s12]
            s33 = nu * (s(1) + s(2));
            e12_tensor = 0.5 * eeng(3);

            wdv = elemData(e).w(g) * elemData(e).detJ0(g) * THICKNESS;

            Sacc = Sacc + wdv * [s(1); s(2); s33; s(3)];
            Eacc = Eacc + wdv * [eeng(1); eeng(2); 0.0; e12_tensor];

            ALLSE = ALLSE + wdv * 0.5 * (s.' * eeng);

            Xr = elemData(e).Xgp(g) - Xc;
            Yr = elemData(e).Ygp(g) - Yc;
            svec = [s(1); s(2); s(3)];     % S11, S22, S12
            Q1acc = Q1acc + wdv * svec * Xr;
            Q2acc = Q2acc + wdv * svec * Yr;
        end
    end

    Savg = Sacc / V0;
    Eavg = Eacc / V0;
    Q1_S = Q1acc / V0;
    Q2_S = Q2acc / V0;
end


% ============================================================
% DIAGNOSTICS -- same Gauss-quadrature / boundary-reaction
% cross-check philosophy as RUC_FE_2D.m's postprocessFrame (Hill
% residual, Q-symmetry residual, global force balance, KUBC recovery
% error), ported to linear elasticity. Two checks have no finite-strain
% analogue and are replaced: EnergyStressResid (which cross-checked a
% hand-coded stress against a finite-difference energy derivative --
% meaningless here since S=D*eps and W=0.5*eps'*D*eps are trivially
% consistent by construction) becomes EnergyConsistencyResid, the
% standard linear-FE identity ALLSE == 0.5*u'*K*u -- a genuine assembly
% self-check, since ALLSE (Gauss quadrature) and K (also Gauss
% quadrature, same B/D matrices) are built independently in this code.
% Jmin (deformation-gradient positivity) has no meaning at all in small
% strain and is simply omitted.
% ============================================================

function [HillResid, QSymResid, ForceBalanceResid, KubcErrF, KubcErrG, EnergyConsistencyResid] = ...
        computeDiagnostics_LE(Rfull, nodes, bndAll, Xc, Yc, Savg, Q1_S, Q2_S, V0, ...
        Fbar, Gm, H, G, u, K, ALLSE)

    Q1acc = Q1_S * V0;
    Q2acc = Q2_S * V0;

    HillPV0 = zeros(2,2);
    QB_xxx=0; QB_xyy=0; QB_xxy=0;
    QB_yxx=0; QB_yyy=0; QB_yxy=0;
    sumRFx=0; sumRFy=0;
    for idx = 1:numel(bndAll)
        nlab = bndAll(idx);
        xr = nodes(nlab,1)-Xc; yr = nodes(nlab,2)-Yc;
        Fx = Rfull(2*nlab-1); Fy = Rfull(2*nlab);
        HillPV0(1,1)=HillPV0(1,1)+Fx*xr; HillPV0(1,2)=HillPV0(1,2)+Fx*yr;
        HillPV0(2,1)=HillPV0(2,1)+Fy*xr; HillPV0(2,2)=HillPV0(2,2)+Fy*yr;
        QB_xxx=QB_xxx+Fx*xr*xr; QB_xyy=QB_xyy+Fx*yr*yr; QB_xxy=QB_xxy+Fx*xr*yr;
        QB_yxx=QB_yxx+Fy*xr*xr; QB_yyy=QB_yyy+Fy*yr*yr; QB_yxy=QB_yxy+Fy*xr*yr;
        sumRFx=sumRFx+Fx; sumRFy=sumRFy+Fy;
    end

    % ---- Hill / macro-homogeneity residual: domain-Gauss <S> (times
    %      V0) versus the boundary-reaction first moment ----
    Snum2x2 = V0 * [Savg(1), Savg(4); Savg(4), Savg(2)];
    hillScale = max(1e-30, max(abs(Snum2x2(:))));
    HillResid = max(abs(Snum2x2(:) - HillPV0(:))) / hillScale;

    % ---- double-stress symmetry residual: domain-Gauss Q (linear-in-X
    %      moment of S) versus the boundary-reaction second moment ----
    domSym_1_11 = 2*Q1acc(1);
    domSym_1_22 = 2*Q2acc(3);
    domSym_1_12 = Q2acc(1) + Q1acc(3);
    domSym_2_11 = 2*Q1acc(3);
    domSym_2_22 = 2*Q2acc(2);
    domSym_2_12 = Q2acc(3) + Q1acc(2);
    domSymVals = [domSym_1_11, domSym_1_22, domSym_1_12, domSym_2_11, domSym_2_22, domSym_2_12];
    QBVals     = [QB_xxx, QB_xyy, QB_xxy, QB_yxx, QB_yyy, QB_yxy];
    qScale = max(1e-30, max(abs(domSymVals)));
    QSymResid = max(abs(domSymVals - QBVals)) / qScale;

    % ---- global force balance: with no body force, the DEDUPLICATED
    %      sum of nodal reactions over every boundary node must vanish ----
    forceScale = max(1e-30, max(abs(Rfull([2*bndAll-1; 2*bndAll]))));
    ForceBalanceResid = norm([sumRFx, sumRFy]) / forceScale;

    % ---- KUBC recovery error: Method-A Fbar/G versus the H,G that was
    %      actually imposed as the Dirichlet data ----
    KubcErrF = max(max(abs(Fbar - (eye(2) + H))));
    KubcErrG = max(abs(Gm(:) - G(:)));

    % ---- energy-assembly consistency: ALLSE (Gauss quadrature) versus
    %      the standard linear-FE identity 0.5*u'*K*u ----
    ALLSE_check = 0.5 * (u.' * (K*u));
    EnergyConsistencyResid = abs(ALLSE - ALLSE_check) / max(1e-30, abs(ALLSE));
end


% ============================================================
% Uavg: area(volume)-weighted mean-of-6-nodal-U per element
% (matches Main2D_LE.py's corrected, non-node-count Uavg)
% ============================================================

function Uavg = elementAveragedU(u, elems6, elemData, THICKNESS)
    Nelem = size(elems6,1);
    u1acc = 0.0; u2acc = 0.0; Vacc = 0.0;
    for e = 1:Nelem
        nodesE = elems6(e,:);
        eVol = sum(elemData(e).w .* elemData(e).detJ0) * THICKNESS;
        u1sum = 0.0; u2sum = 0.0;
        for k = 1:6
            n = nodesE(k);
            u1sum = u1sum + u(2*n-1);
            u2sum = u2sum + u(2*n);
        end
        u1acc = u1acc + eVol * (u1sum/6);
        u2acc = u2acc + eVol * (u2sum/6);
        Vacc  = Vacc + eVol;
    end
    Uavg = [u1acc; u2acc] / Vacc;
end


% ============================================================
% NODAL (SMOOTHED) FIELDS FOR PLOTTING  -- Sij, Eij, Ui evaluated AT
% each element's own 6 nodal parametric positions (not just the 3
% Gauss points used for stiffness/domain integration), then arithmetic-
% averaged over every element sharing that global node. Standard FE
% "smoothed stress" post-processing -- needed because Sij/Eij are only
% C(-1) (discontinuous) across element edges otherwise, which plots
% poorly. U1,U2 are primary nodal DOFs already, so no averaging is
% needed for them -- taken directly from u.
% ============================================================

function [Sfield, Efield, Ufield] = computeNodalFields(u, nodes, elems6, Delem, nu_e)
    Nnodes = size(nodes,1);
    Nelem  = size(elems6,1);

    % natural (L1,L2) coordinates of the 6 T6 nodes, Abaqus/PDE-Toolbox
    % order: 3 corners, then 3 midsides on edges 1-2, 2-3, 3-1.
    nodeL1 = [1, 0, 0, 0.5, 0.0, 0.5];
    nodeL2 = [0, 1, 0, 0.5, 0.5, 0.0];

    Sacc = zeros(Nnodes,4);   % S11,S22,S33,S12
    Eacc = zeros(Nnodes,4);   % E11,E22,E33(=0),E12
    cnt  = zeros(Nnodes,1);

    for e = 1:Nelem
        nodesE = elems6(e,:);
        coords6 = nodes(nodesE,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        Dmat = Delem{e};
        nu = nu_e(e);

        for k = 1:6
            [~, dNdX] = t6ShapeDeriv(nodeL1(k), nodeL2(k), coords6);
            B = zeros(3,12);
            for a = 1:6
                B(1, 2*a-1) = dNdX(a,1);
                B(2, 2*a)   = dNdX(a,2);
                B(3, 2*a-1) = dNdX(a,2);
                B(3, 2*a)   = dNdX(a,1);
            end
            eeng = B * ue;
            s = Dmat * eeng;
            s33 = nu * (s(1) + s(2));
            e12 = 0.5 * eeng(3);

            n = nodesE(k);
            Sacc(n,:) = Sacc(n,:) + [s(1), s(2), s33, s(3)];
            Eacc(n,:) = Eacc(n,:) + [eeng(1), eeng(2), 0.0, e12];
            cnt(n) = cnt(n) + 1;
        end
    end

    Sfield = Sacc ./ cnt;
    Efield = Eacc ./ cnt;
    Ufield = [u(1:2:end), u(2:2:end)];   % Nnodes x 2, primary DOFs, no averaging
end


% ============================================================
% FIELD PLOTS: Sij, Eij, Ui on the mesh ("plot function with msh")
% Uses the corner-triangle sub-connectivity of each T6 element for
% `patch` (same simplification RUC_FE_2D.m's own quick-look plot
% uses); colour is the smoothed nodal field, interpolated across each
% triangle. Plotted on the UNDEFORMED mesh by default (deformScale=0);
% pass deformScale>0 to exaggerate the (tiny, ~strain0) displacement.
% ============================================================

function plotStressStrainDisp(nodes, elems6, Sfield, Efield, Ufield, deformScale, tagStr)

    tri = elems6(:,1:3);
    Xp = nodes(:,1) + deformScale*Ufield(:,1);
    Yp = nodes(:,2) + deformScale*Ufield(:,2);

    figure('Color','w', 'Name', ['Sij  --  ' tagStr], 'Position',[80 80 1000 800]);
    Stitles = {'S_{11}','S_{22}','S_{33}','S_{12}'};
    for k = 1:4
        subplot(2,2,k);
        patch('Faces',tri, 'Vertices',[Xp,Yp], 'FaceVertexCData',Sfield(:,k), ...
              'FaceColor','interp', 'EdgeColor','none');
        axis equal tight; colorbar; title(Stitles{k});
    end
    sgtitle(['Stress field -- ' tagStr], 'Interpreter','none');

    figure('Color','w', 'Name', ['Eij  --  ' tagStr], 'Position',[120 60 1000 800]);
    Etitles = {'E_{11}','E_{22}','E_{33}','E_{12}'};
    for k = 1:4
        subplot(2,2,k);
        patch('Faces',tri, 'Vertices',[Xp,Yp], 'FaceVertexCData',Efield(:,k), ...
              'FaceColor','interp', 'EdgeColor','none');
        axis equal tight; colorbar; title(Etitles{k});
    end
    sgtitle(['Strain field -- ' tagStr], 'Interpreter','none');

    figure('Color','w', 'Name', ['U  --  ' tagStr], 'Position',[160 40 1000 460]);
    Utitles = {'U_{1}','U_{2}'};
    for k = 1:2
        subplot(1,2,k);
        patch('Faces',tri, 'Vertices',[Xp,Yp], 'FaceVertexCData',Ufield(:,k), ...
              'FaceColor','interp', 'EdgeColor','none');
        axis equal tight; colorbar; title(Utitles{k});
    end
    sgtitle(['Displacement field -- ' tagStr], 'Interpreter','none');
end


% ============================================================
% REACTION-FORCE CURVE PLOTS on LEFT/RIGHT/BOTTOM/TOP -- the actual
% nodal Rfull distribution along each side (not just its RF_*_x/y sum),
% plotted against position along that side.
% ============================================================

function plotReactionForceCurves(Rfull, nodes, leftN, rightN, botN, topN, tagStr)

    [yL, FxL, FyL] = sideCurve(Rfull, nodes, leftN,  2);
    [yR, FxR, FyR] = sideCurve(Rfull, nodes, rightN, 2);
    [xB, FxB, FyB] = sideCurve(Rfull, nodes, botN,   1);
    [xT, FxT, FyT] = sideCurve(Rfull, nodes, topN,   1);

    figure('Color','w', 'Name', ['Reaction forces  --  ' tagStr], 'Position',[200 100 1000 700]);

    subplot(2,2,1);
    plot(yL,FxL,'-o', yR,FxR,'-s', 'LineWidth',1.4, 'MarkerSize',5);
    grid on; xlabel('Y'); ylabel('R_x'); title('LEFT / RIGHT: R_x(Y)');
    legend({'LEFT','RIGHT'}, 'Location','best');

    subplot(2,2,2);
    plot(yL,FyL,'-o', yR,FyR,'-s', 'LineWidth',1.4, 'MarkerSize',5);
    grid on; xlabel('Y'); ylabel('R_y'); title('LEFT / RIGHT: R_y(Y)');
    legend({'LEFT','RIGHT'}, 'Location','best');

    subplot(2,2,3);
    plot(xB,FxB,'-o', xT,FxT,'-s', 'LineWidth',1.4, 'MarkerSize',5);
    grid on; xlabel('X'); ylabel('R_x'); title('BOTTOM / TOP: R_x(X)');
    legend({'BOTTOM','TOP'}, 'Location','best');

    subplot(2,2,4);
    plot(xB,FyB,'-o', xT,FyT,'-s', 'LineWidth',1.4, 'MarkerSize',5);
    grid on; xlabel('X'); ylabel('R_y'); title('BOTTOM / TOP: R_y(X)');
    legend({'BOTTOM','TOP'}, 'Location','best');

    sgtitle(['Boundary reaction-force distribution -- ' tagStr], 'Interpreter','none');
end

function [pos, Fx, Fy] = sideCurve(Rfull, nodes, nodeList, sortCol)
    pos = nodes(nodeList, sortCol);
    Fx  = Rfull(2*nodeList-1);
    Fy  = Rfull(2*nodeList);
    [pos, ord] = sort(pos);
    Fx = Fx(ord); Fy = Fy(ord);
end


% ============================================================
% Fbar/G RECOVERY -- METHOD A: boundary least-squares fit
% (material-independent, purely geometric -- same basis as
% Main2D_LE.py's Method A and RUC_FE_2D.m's solveGauss/fit)
% ============================================================

function [Fbar, G, fitResidMax, c1_fit, c2_fit] = boundaryFbarG_LSQ( ...
        u, nodes, bndAll, Xc, Yc)

    nB = numel(bndAll);
    Phi = zeros(nB,6);
    U1 = zeros(nB,1); U2 = zeros(nB,1);
    for k = 1:nB
        n = bndAll(k);
        p = nodes(n,1) - Xc;
        q = nodes(n,2) - Yc;
        Phi(k,:) = [1, p, q, 0.5*p*p, p*q, 0.5*q*q];
        U1(k) = u(2*n-1);
        U2(k) = u(2*n);
    end

    coef1 = Phi \ U1;   % [c1, H11, H12, G111, G112(=G121), G122]
    coef2 = Phi \ U2;   % [c2, H21, H22, G211, G212(=G221), G222]

    c1_fit = coef1(1); c2_fit = coef2(1);

    Fbar = [1 + coef1(2), coef1(3); ...
            coef2(2),     1 + coef2(3)];

    G = zeros(2,2,2);
    G(1,1,1) = coef1(4); G(1,1,2) = coef1(5); G(1,2,1) = coef1(5); G(1,2,2) = coef1(6);
    G(2,1,1) = coef2(4); G(2,1,2) = coef2(5); G(2,2,1) = coef2(5); G(2,2,2) = coef2(6);

    pred1 = Phi * coef1; pred2 = Phi * coef2;
    fitResidMax = max( max(abs(pred1 - U1)), max(abs(pred2 - U2)) );
end


% ============================================================
% Fbar/G RECOVERY -- METHOD B: corrected surface integral
% (trapezoidal edge quadrature; same closed-form correction as
% Main2D_LE.py / RUC_FE_2D.m's boundaryFbarG_FE_matlab)
% ============================================================

function w = trapWeights1D(p)
    n = numel(p);
    w = zeros(n,1);
    if n == 1, return; end
    w(1) = 0.5*(p(2)-p(1));
    w(n) = 0.5*(p(n)-p(n-1));
    for i = 2:n-1
        w(i) = 0.5*(p(i+1)-p(i-1));
    end
end

function [dat, w] = edgeData(u, nodes, nlist, sortCol)
    x = nodes(nlist,1); y = nodes(nlist,2);
    u1 = u(2*nlist-1); u2 = u(2*nlist);
    dat = [x, y, u1, u2];
    [~, ord] = sort(dat(:,sortCol));
    dat = dat(ord,:);
    w = trapWeights1D(dat(:,sortCol));
end

function [Fbar, G] = boundaryFbarG_Integral(u, nodes, leftN, rightN, botN, topN, ...
        Lx_tot, Ly_tot, Xc, Yc)

    Vtot = Lx_tot * Ly_tot;
    [left, wL]   = edgeData(u, nodes, leftN,  2);
    [right, wR]  = edgeData(u, nodes, rightN, 2);
    [bottom, wB] = edgeData(u, nodes, botN,   1);
    [top, wT]    = edgeData(u, nodes, topN,   1);

    u1L=left(:,3); u2L=left(:,4); ycL=left(:,2)-Yc;
    u1R=right(:,3);u2R=right(:,4);ycR=right(:,2)-Yc;
    u1B=bottom(:,3);u2B=bottom(:,4);xcB=bottom(:,1)-Xc;
    u1T=top(:,3);   u2T=top(:,4);   xcT=top(:,1)-Xc;

    L_Iu1=wL'*u1L; L_Iu2=wL'*u2L;
    R_Iu1=wR'*u1R; R_Iu2=wR'*u2R;
    B_Iu1=wB'*u1B; B_Iu2=wB'*u2B;
    T_Iu1=wT'*u1T; T_Iu2=wT'*u2T;

    xcL_c = left(1,1)-Xc;  xcR_c = right(1,1)-Xc;
    ycB_c = bottom(1,2)-Yc; ycT_c = top(1,2)-Yc;

    L_Iu1x1=xcL_c*L_Iu1; L_Iu1x2=wL'*(u1L.*ycL);
    L_Iu2x1=xcL_c*L_Iu2; L_Iu2x2=wL'*(u2L.*ycL);
    R_Iu1x1=xcR_c*R_Iu1; R_Iu1x2=wR'*(u1R.*ycR);
    R_Iu2x1=xcR_c*R_Iu2; R_Iu2x2=wR'*(u2R.*ycR);
    B_Iu1x1=wB'*(u1B.*xcB); B_Iu1x2=ycB_c*B_Iu1;
    B_Iu2x1=wB'*(u2B.*xcB); B_Iu2x2=ycB_c*B_Iu2;
    T_Iu1x1=wT'*(u1T.*xcT); T_Iu1x2=ycT_c*T_Iu1;
    T_Iu2x1=wT'*(u2T.*xcT); T_Iu2x2=ycT_c*T_Iu2;

    F11 = 1 + (-L_Iu1+R_Iu1)/Vtot;
    F12 = (-B_Iu1+T_Iu1)/Vtot;
    F21 = (-L_Iu2+R_Iu2)/Vtot;
    F22 = 1 + (-B_Iu2+T_Iu2)/Vtot;

    S111=(-L_Iu1x1+R_Iu1x1)/Vtot; S112=(-L_Iu1x2+R_Iu1x2)/Vtot;
    S121=(-B_Iu1x1+T_Iu1x1)/Vtot; S122=(-B_Iu1x2+T_Iu1x2)/Vtot;
    S211=(-L_Iu2x1+R_Iu2x1)/Vtot; S212=(-L_Iu2x2+R_Iu2x2)/Vtot;
    S221=(-B_Iu2x1+T_Iu2x1)/Vtot; S222=(-B_Iu2x2+T_Iu2x2)/Vtot;

    varX = Lx_tot^2/12; varY = Ly_tot^2/12;
    % closed-form mean of the PRESCRIBED polynomial -- boundary data
    % only, NOT the FE field's own <u>_V (see project learnings: this
    % is the volume-averaging correction that fixed G111/G122/G211/G222).
    A1 = 0.25*(S111+S122); A2 = 0.25*(S211+S222);

    G111=(S111-A1)/varX; G112=S112/varY;
    G121=S121/varX;      G122=(S122-A1)/varY;
    G211=(S211-A2)/varX; G212=S212/varY;
    G221=S221/varX;      G222=(S222-A2)/varY;

    Fbar = [F11 F12; F21 F22];
    G = zeros(2,2,2);
    G(1,1,1)=G111; G(1,1,2)=G112; G(1,2,1)=G121; G(1,2,2)=G122;
    G(2,1,1)=G211; G(2,1,2)=G212; G(2,2,1)=G221; G(2,2,2)=G222;
end


% ============================================================
% PER-SIDE BOUNDARY REACTION FORCE (Rfull already includes zero at
% every free dof by construction of the linear solve)
% ============================================================

function [Fx, Fy] = sideReactionForces(Rfull, nodeList)
    Fx = sum(Rfull(2*nodeList-1));
    Fy = sum(Rfull(2*nodeList));
end


% ============================================================
% GEOMETRY / MESH  (reused unchanged from RUC_FE_2D.m -- purely
% geometric, material-independent)
% ============================================================

function tf = circleIntersectsDomain(cx,cy,Rc,xmin,xmax,ymin,ymax)
    tf = ~(cx+Rc < xmin || cx-Rc > xmax || cy+Rc < ymin || cy-Rc > ymax);
end

function centers = buildInclusionCenters(NX,NY,Lx,Ly,zeta1,zeta2,Lx_tot,Ly_tot,R)
    centers = zeros(0,2);
    for i = -1:NX
        for j = -1:NY
            cx = (i+0.5)*Lx + zeta1;
            cy = (j+0.5)*Ly + zeta2;
            if circleIntersectsDomain(cx,cy,R,0,Lx_tot,0,Ly_tot)
                centers(end+1,:) = [cx,cy]; %#ok<AGROW>
            end
        end
    end
    if isempty(centers)
        error('RUC_FE_2D_LinearElastic:NoInclusions', ...
            'No inclusion circles intersect the domain -- check NX/NY/zeta/geometry.');
    end
end

function [nodes, elems6, isInclusionElem] = buildRUCMesh(Lx_tot,Ly_tot,centers,R,meshSize)
    if exist('createpde','file') ~= 2
        error('RUC_FE_2D_LinearElastic:NoPDEToolbox', ['This mesh generator needs the ' ...
            'Partial Differential Equation Toolbox (createpde/geometryFromEdges/' ...
            'generateMesh). If unavailable, supply your own T6 mesh: set nodes ' ...
            '(Nx2) and elems6 (Mx6, 3 corners then 3 midsides per Abaqus/PDE-' ...
            'Toolbox node order) directly.']);
    end
    Nc = size(centers,1);
    rectGd = [3;4; 0;Lx_tot;Lx_tot;0; 0;0;Ly_tot;Ly_tot];
    gd = rectGd;
    names = {'R1'};
    sf = 'R1';
    for k = 1:Nc
        circGd = [1; centers(k,1); centers(k,2); R; 0;0;0;0;0;0];
        gd = [gd, circGd]; %#ok<AGROW>
        nm = sprintf('C%d',k);
        names{end+1} = nm; %#ok<AGROW>
        sf = [sf '+' nm]; %#ok<AGROW>
    end
    ns = char(names)';
    dl = decsg(gd, sf, ns);

    model = createpde(1);
    geometryFromEdges(model, dl);
    generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'quadratic');

    nodes = model.Mesh.Nodes';
    elems6 = model.Mesh.Elements';

    Nelem = size(elems6,1);
    isInclusionElem = false(Nelem,1);
    for e = 1:Nelem
        cn = elems6(e,1:3);
        cx = mean(nodes(cn,1)); cy = mean(nodes(cn,2));
        d2 = (centers(:,1)-cx).^2 + (centers(:,2)-cy).^2;
        isInclusionElem(e) = any(d2 <= R^2);
    end
end

function [leftN,rightN,botN,topN,bndAll] = classifyBoundary(nodes,Lx,Ly,TOL)
    x = nodes(:,1); y = nodes(:,2);
    leftN  = find(abs(x-0)  < TOL);
    rightN = find(abs(x-Lx) < TOL);
    botN   = find(abs(y-0)  < TOL);
    topN   = find(abs(y-Ly) < TOL);
    bndAll = unique([leftN;rightN;botN;topN]);
end

function [N, dNdX, detJ, X, Y] = t6ShapeDeriv(L1,L2,coords6)
    L3 = 1-L1-L2;
    N = [L1*(2*L1-1); L2*(2*L2-1); L3*(2*L3-1); 4*L1*L2; 4*L2*L3; 4*L3*L1];
    dN_dL1 = [4*L1-1; 0; -(4*L3-1); 4*L2; -4*L2; 4*(L3-L1)];
    dN_dL2 = [0; 4*L2-1; -(4*L3-1); 4*L1; 4*(L3-L2); -4*L1];
    x = coords6(:,1); y = coords6(:,2);
    dXdL1 = dN_dL1'*x; dXdL2 = dN_dL2'*x;
    dYdL1 = dN_dL1'*y; dYdL2 = dN_dL2'*y;
    J = [dXdL1 dXdL2; dYdL1 dYdL2];
    detJ = det(J);
    Jinv = inv(J);
    dL1dX = Jinv(1,1); dL2dX = Jinv(2,1);
    dL1dY = Jinv(1,2); dL2dY = Jinv(2,2);
    dNdX = zeros(6,2);
    dNdX(:,1) = dN_dL1*dL1dX + dN_dL2*dL2dX;
    dNdX(:,2) = dN_dL1*dL1dY + dN_dL2*dL2dY;
    X = N'*x; Y = N'*y;
end

function [L1w,L2w,W] = gauss3()
    L1w = [2/3, 1/6, 1/6];
    L2w = [1/6, 2/3, 1/6];
    W   = [1/6, 1/6, 1/6];
end


% ============================================================
% HOMOGENEOUS (C_hom) REFERENCE SOLVE -- NEW
% Plain single-material mesh (no inclusion) and its per-load-case
% KUBC solve/domain-average, used to fill the *_hom CSV columns.
% ============================================================

function [nodes, elems6] = buildPlainMesh(Lx_tot, Ly_tot, meshSize)
    if exist('createpde','file') ~= 2
        error('RUC_FE_2D_LinearElastic:NoPDEToolbox', ['buildPlainMesh needs the ' ...
            'Partial Differential Equation Toolbox (createpde/geometryFromEdges/' ...
            'generateMesh).']);
    end
    gd = [3;4; 0;Lx_tot;Lx_tot;0; 0;0;Ly_tot;Ly_tot];
    ns = char({'R1'})';
    sf = 'R1';
    dl = decsg(gd, sf, ns);

    model = createpde(1);
    geometryFromEdges(model, dl);
    generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'quadratic');

    nodes = model.Mesh.Nodes';
    elems6 = model.Mesh.Elements';
end

function [Savg, Eavg, ALLSE, Q1_S, Q2_S] = domainAverages_Hom(u, elems6, elemData, Dmat, Cs33row, THICKNESS, V0, Xc, Yc)
    % Same [11,22,33,12] slot convention and Q-moment definition as
    % domainAverages_LE, but for a single homogeneous Dmat: S33 has no
    % per-element "nu" to come from, so it's filled in AFTER the loop
    % from the homogenized out-of-plane row Cs33row. Energy is
    % unaffected by that choice: in plane strain E33=0, so S33*E33
    % never contributes to ALLSE regardless of S33's value.
    Nelem = size(elems6,1);
    Sacc = zeros(4,1); Eacc = zeros(4,1); ALLSE = 0.0;
    Q1acc = zeros(3,1); Q2acc = zeros(3,1);
    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        for g = 1:3
            dNdX = elemData(e).dNdX{g};
            B = zeros(3,12);
            for a = 1:6
                B(1, 2*a-1) = dNdX(a,1);
                B(2, 2*a)   = dNdX(a,2);
                B(3, 2*a-1) = dNdX(a,2);
                B(3, 2*a)   = dNdX(a,1);
            end
            eeng = B * ue;               % [e11; e22; gamma12]
            s = Dmat * eeng;             % [s11; s22; s12]
            wdv = elemData(e).w(g) * elemData(e).detJ0(g) * THICKNESS;

            Sacc = Sacc + wdv * [s(1); s(2); 0; s(3)];
            Eacc = Eacc + wdv * [eeng(1); eeng(2); 0; 0.5*eeng(3)];
            ALLSE = ALLSE + wdv * 0.5 * (s.' * eeng);

            Xr = elemData(e).Xgp(g) - Xc;
            Yr = elemData(e).Ygp(g) - Yc;
            svec = [s(1); s(2); s(3)];
            Q1acc = Q1acc + wdv * svec * Xr;
            Q2acc = Q2acc + wdv * svec * Yr;
        end
    end
    Savg = Sacc / V0;
    Eavg = Eacc / V0;
    Savg(3) = Cs33row * [Eavg(1); Eavg(2); 2*Eavg(4)];   % fill S33 post-hoc
    Q1_S = Q1acc / V0;
    Q2_S = Q2acc / V0;
end

function row = solveHomogeneousCase( ...
        K_hom, nodesH, elems6H, elemDataH, leftNH,rightNH,botNH,topNH,bndAllH, ...
        C_hom, Cs33_hom, THICKNESS, V0_hom, Lx_tot, Ly_tot, Xc, Yc, ndofH, ...
        H, G, caseName, strain0, caseIdx)
    % One full KUBC solve + the SAME post-processing as the
    % heterogeneous per-case block, but on the plain homogeneous mesh
    % with material C_hom -- returns a row struct with exactly the
    % heterogeneous table's columns (minus NX/NY/zeta1/zeta2).

    [isDirichlet, uBC] = buildKUBC_BC(H, G, nodesH, bndAllH, Xc, Yc, ndofH);
    freeDofs = find(~isDirichlet);
    dirDofs  = find(isDirichlet);
    u = zeros(ndofH,1);
    u(dirDofs) = uBC(dirDofs);
    if ~isempty(freeDofs)
        Kff = K_hom(freeDofs, freeDofs);
        Kfd = K_hom(freeDofs, dirDofs);
        rhs = -Kfd * u(dirDofs);
        u(freeDofs) = Kff \ rhs;
    end
    Rfull = K_hom * u;

    [RF_L_x, RF_L_y] = sideReactionForces(Rfull, leftNH);
    [RF_R_x, RF_R_y] = sideReactionForces(Rfull, rightNH);
    [RF_B_x, RF_B_y] = sideReactionForces(Rfull, botNH);
    [RF_T_x, RF_T_y] = sideReactionForces(Rfull, topNH);

    [Savg, Eavg, ALLSE, Q1_S, Q2_S] = domainAverages_Hom( ...
        u, elems6H, elemDataH, C_hom, Cs33_hom, THICKNESS, V0_hom, Xc, Yc);

    Q111 = Q1_S(1); Q112 = Q2_S(1);
    Q221 = Q1_S(2); Q222 = Q2_S(2);
    Q121 = Q1_S(3); Q122 = Q2_S(3);
    Q211 = Q1_S(3); Q212 = Q2_S(3);

    Uavg = elementAveragedU(u, elems6H, elemDataH, THICKNESS);

    [Fbar, Gm, fitResidMax, c1_fit, c2_fit] = boundaryFbarG_LSQ( ...
        u, nodesH, bndAllH, Xc, Yc);
    [Fbar_int, G_int] = boundaryFbarG_Integral( ...
        u, nodesH, leftNH, rightNH, botNH, topNH, Lx_tot, Ly_tot, Xc, Yc);

    [HillResid, QSymResid, ForceBalanceResid, KubcErrF, KubcErrG, EnergyConsistencyResid] = ...
        computeDiagnostics_LE(Rfull, nodesH, bndAllH, Xc, Yc, Savg, Q1_S, Q2_S, V0_hom, ...
            Fbar, Gm, H, G, u, K_hom, ALLSE);

    row = struct();
    row.Case = caseIdx;
    row.LoadCase = caseName;
    row.strain0 = strain0;
    row.StepTime = 1.0;
    row.ALLSE = ALLSE;
    row.RF_LEFT_x = RF_L_x;   row.RF_LEFT_y = RF_L_y;
    row.RF_RIGHT_x = RF_R_x;  row.RF_RIGHT_y = RF_R_y;
    row.RF_BOTTOM_x = RF_B_x; row.RF_BOTTOM_y = RF_B_y;
    row.RF_TOP_x = RF_T_x;    row.RF_TOP_y = RF_T_y;
    row.Volume = V0_hom;
    row.nBndNodes = numel(bndAllH);
    row.Savg_11 = Savg(1); row.Savg_22 = Savg(2);
    row.Savg_33 = Savg(3); row.Savg_12 = Savg(4);
    row.Eavg_11 = Eavg(1); row.Eavg_22 = Eavg(2);
    row.Eavg_33 = Eavg(3); row.Eavg_12 = Eavg(4);
    row.Uavg_1 = Uavg(1);  row.Uavg_2 = Uavg(2);
    row.Fbar_11 = Fbar(1,1); row.Fbar_12 = Fbar(1,2);
    row.Fbar_21 = Fbar(2,1); row.Fbar_22 = Fbar(2,2);
    row.G111 = Gm(1,1,1); row.G112 = Gm(1,1,2);
    row.G121 = Gm(1,2,1); row.G122 = Gm(1,2,2);
    row.G211 = Gm(2,1,1); row.G212 = Gm(2,1,2);
    row.G221 = Gm(2,2,1); row.G222 = Gm(2,2,2);
    row.Fbar_11_int = Fbar_int(1,1); row.Fbar_12_int = Fbar_int(1,2);
    row.Fbar_21_int = Fbar_int(2,1); row.Fbar_22_int = Fbar_int(2,2);
    row.G111_int = G_int(1,1,1); row.G112_int = G_int(1,1,2);
    row.G121_int = G_int(1,2,1); row.G122_int = G_int(1,2,2);
    row.G211_int = G_int(2,1,1); row.G212_int = G_int(2,1,2);
    row.G221_int = G_int(2,2,1); row.G222_int = G_int(2,2,2);
    row.Q111 = Q111; row.Q112 = Q112; row.Q121 = Q121; row.Q122 = Q122;
    row.Q211 = Q211; row.Q212 = Q212; row.Q221 = Q221; row.Q222 = Q222;
    row.fitResidMax = fitResidMax; row.c1_fit = c1_fit; row.c2_fit = c2_fit;
    row.ALLSE_std = 0.0;   % no translation ensemble for the homogeneous solve
    row.HillResid = HillResid; row.QSymResid = QSymResid;
    row.ForceBalanceResid = ForceBalanceResid;
    row.KubcErrF = KubcErrF; row.KubcErrG = KubcErrG;
    row.EnergyConsistencyResid = EnergyConsistencyResid;
end