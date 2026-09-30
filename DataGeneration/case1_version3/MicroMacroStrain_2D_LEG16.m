% ============================================================
%  MicroMacroStrain_2D_LE.m
%
%  Macro strain eps_M and translation-averaged micro strain sym(psi)
%  from classical 2D linear-elastic FE solves, following the note
%  "Extraction of the Micro-Strain Field sym(psi) from Finite-Element
%  Solves -- A Three-Stage Procedure with Translation Averaging":
%
%   STAGE 1  plain domain [0,Lx_tot]x[0,Ly_tot], material C_hom from
%            homogenize2D_PBC.m, KUBC on every boundary node
%               u_i = H_ij Xc_j + 1/2 G_ijk Xc_j Xc_k ,  Xc = X - X_c
%            -> macro strain eps_XX, eps_YY, gamma_XY over the domain
%
%   STAGE 2  heterogeneous N x N RUC (matrix + circular inclusions),
%            SAME KUBC, inclusions translated by zeta^(j), j=1..Nzeta.
%            For sym(psi) -> eps_M as N grows, the shifts must cover a
%            whole unit cell (translationMode = 5); shifts over a part of
%            the cell leave the periodic micro-fluctuation in sym(psi).
%            The strain of every solve is evaluated on ONE fixed grid
%            (the meshes differ between translations) and averaged over
%            the translations, separately for every load case k:
%               sym(psi)^(k)(X) = 1/Nzeta * sum_j eps^(k,j)(X)
%
%   STAGE 3  plots of sym(psi)^(k) and of
%               Delta^(k) = sym(psi)^(k) - eps_M
%
%   STAGE 4  relative strain and curvature
%               gamma = eps_M - sym(psi)          (= -Delta)
%               kappa = grad sym(psi), 6-vector
%                 [k_xxx k_xxy k_yyx k_yyy k_xyx k_xyy]
%               = [d psi11/dX, d psi11/dY, d psi22/dX, d psi22/dY,
%                  d psi12/dX, d psi12/dY]      (tensor psi12)
%            and the domain integrals of the energy terms, per load case
%               T1 = int eps_M' C eps_M      T2 = int sym(psi)' C sym(psi)
%               T3 = int gamma' C gamma      T4 = int kappa' A_kappa kappa
%               T5 = int sym(psi)' C gamma   (coupling term of W)
%            with C = C_hom (Voigt, engineering shear), so that
%               W = am/2 T2 + ah/2 T3 + ac T5 + al/2 (L/N)^2 T4 ,
%            and T1 = T2 + 2 T5 + T3 identically (gamma = eps_M - psi).
%            ORTHOGONALITY (macroRef = 'stage1FE'): T3 + T5 = int eps_M' C gamma
%            = 0 in the continuum, because C eps_M is equilibrated (Stage 1)
%            and gamma = sym grad(u_M - mean_j u^(j)) with u_M - mean_j u^(j)
%            = 0 on the whole KUBC boundary. Hence T5 = -T3, T2 = T1 + T3 and
%               W = am/2 T1 + (am + ah - 2 ac)/2 T3 + al/2 (L/N)^2 T4 ,
%            i.e. the energy sees alpha_h, alpha_c only through
%            am + ah - 2 ac. The script prints the discrete residual.
%            Reference energies: 2*ALLSE of the Stage-1 solve (= T1 up to
%            grid quadrature) and 2*ALLSE of the heterogeneous solves.
%
%  Conventions are identical to RUC_FE_2D_LinearElastic.m: raw 9-vector
%  [H11,H22,H12,G1_11,G1_22,G1_12,G2_11,G2_22,G2_12], H = raw(1:3)*strain0,
%  G = raw(4:9)*strain0*gScaleFactor, CPE6-type T6 elements with 3-point
%  Gauss rule, plane strain, translationMode 1..4, same mesher.
%
%  NOTE ON THE MACRO REFERENCE (macroRef below). The note states that
%  the Stage-1 field equals sym(H) everywhere. That is exact only when
%  the KUBC polynomial is itself an equilibrium field of C_hom, i.e.
%  div(C_hom : sym(G)) = 0 (always true for G = 0). For a G-load such
%  as raw = [0 0 0 0 0 0 1 0 0] (u_2 = 1/2 G_211 X1^2) it is NOT: the
%  polynomial needs a body force b_2 = -C_hom(3,3)*G_211, so the
%  Stage-1 FE field = polynomial + an equilibrating correction, and
%  sym(H) = 0. The script prints that body force per load case.
%  Default macroRef = 'stage1FE' uses the actual Stage-1 FE field, so
%  Delta vanishes for a homogeneous medium for ANY (H,G), which is the
%  property the note's Stage-3 interpretation relies on.
%
%  REQUIRES: homogenize2D_PBC.m on the MATLAB path; PDE Toolbox
%  (decsg/createpde/generateMesh), used only for mesh generation.
%  Uses local functions in a script (MATLAB R2016b+), sgtitle (R2018b+).
% ============================================================

clear; clc; close all;

%% ============ USER INPUT ============

% ---- specimen (fixed, independent of N) ----
Lx_tot = 1.0;  Ly_tot = 1.0;

% ---- heterogeneous RUC used in Stage 2 ----
Nlist    = [1 2 3 4 5 6 7 8]; % unit cells per side (N x N inclusions); Stages 2-4 run for each N
Rfrac    = 6/19;    % inclusion radius / unit-cell size
meshFrac = 0.05;    % element size / unit-cell size  (Stage 2 meshes)

% ---- Stage-1 plain mesh and PBC homogenization mesh ----
meshFracHom = 0.02; % element size / min(Lx_tot,Ly_tot) (Stage 1) and
                    % element size / cell size in homogenize2D_PBC

% ---- loading amplitude (same convention as RUC_FE_2D_LinearElastic.m) ----
strain0      = 0.01;
gScaleFactor = 2;

% ---- material, plane strain ----
E_m = 70000.0;  nu_m = 0.33;   % matrix
E_i = 3500.0;   nu_i = 0.33;   % inclusion  (set E_i=E_m, nu_i=nu_m for the
                               % homogeneous sanity check: Delta -> 0)
THICKNESS = 1.0;

% ------------------------------------------------------------
% LOAD CASE SETTINGS
% ------------------------------------------------------------
loadCaseSource = 'library';  % 'library', 'file' or 'manual'
% used when loadCaseSource = 'library' (raw vectors of unit norm, H = 0):
%   'G6'  : the 6 single G components  G1_11 G1_22 G1_12 G2_11 G2_22 G2_12
%   'G16' : G6 + 10 combinations (see buildCaseLibrary):
%           bending / shear-gradient about X and Y  (G1_12 -/+ G2_11,
%           G2_12 -/+ G1_22), co-/counter-directional normal-strain
%           gradients (G1_11 +/- G2_12, G1_12 +/- G2_22), G1_11+G2_22, all six
%   'G21' : G6 + all 15 pairs (+), i.e. enough to fix every entry of a
%           quadratic form in G
caseLibrary = 'G16';
% used when loadCaseSource = 'file': take the FIRST this many rows of
% sample_points_d09.txt (which must sit next to this .m file). Because
% that file is hierarchically ordered, this IS an optimized N-point set.
numCasesFromFile = 64;
sampleFileName   = 'sample_points_d09.txt';
% used when loadCaseSource = 'manual': one row per case, each a raw
% 9-vector [H11,H22,H12,G1_11,G1_22,G1_12,G2_11,G2_22,G2_12] on the
% UNSCALED (unit-sphere) scale -- multiplied by strain0 below, same as
% the 'file' source. Example: pure H11 stretch direction.
% manualRawVectors = [ ...
%     +5.72245728820545696e-01 -1.92916345731379196e-01 -1.25580798095818208e-01 -4.43331170730187241e-01 -2.71724749117376430e-01 -1.94577047958669619e-01 +5.34194488340966944e-01 -9.46217548967859950e-02 +1.30358930864310052e-01
% ];
manualRawVectors = [0 0 0 0 0 0 1 0 0];

% ------------------------------------------------------------
% TRANSLATION SETTINGS (relative to the CELL; same as the RUC code)
% ------------------------------------------------------------
marginFrac      = 0.05;
stepFrac        = 0.05;
translationMode = 5;   % 1=diagonal, 2=full grid, 3=plus-cross, 4=centred only,
                       % 5=FULL PERIOD: nPeriod x nPeriod uniform shifts spanning
                       %   one whole unit cell (inclusions may cross the outer
                       %   boundary; they are clipped to the domain)
nPeriod         = 8;   % shifts per direction for mode 5 (nPeriod^2 solves per N)
% Modes 1-3 only shift the inclusion by +/-(1/2 - Rfrac - marginFrac) of a
% cell, i.e. over a small part of one period. Their average therefore
% still contains the periodic micro-fluctuation of the strain, whose
% energy does NOT shrink with N, so sym(psi) does not tend to eps_M.
% Mode 5 averages over a whole period, which removes that fluctuation
% while leaving the macroscopic (load-driven) variation untouched.

% ------------------------------------------------------------
% FIXED EVALUATION GRID (common to Stage 1 and every translation)
% ------------------------------------------------------------
nGridX = 121;  nGridY = 121;   % odd -> the midlines X=Xc, Y=Yc are grid lines

% ------------------------------------------------------------
% STAGE-3 MACRO REFERENCE eps_M
%   'stage1FE' : eps_M(X) = strain of the Stage-1 C_hom FE solve (default)
%   'symH'     : eps_M = sym(H), constant (formula written in Stage 3)
%   'affine'   : eps_M(X) = sym(H + G.Xc) (psi ansatz of the identification)
% ------------------------------------------------------------
macroRef = 'stage1FE';

% ------------------------------------------------------------
% STAGE 4: ENERGY TERMS
% ------------------------------------------------------------
% A_kappa = diag(Akappa_diag) acting on
% kappa = [k_xxx k_xxy k_yyx k_yyy k_xyx k_xyy] (same default as
% identify_micromorphic_alphas_vs_N.m). [1 1 1 1 2 2] would instead give
% the full tensor contraction kappa_ijk*kappa_ijk.
Akappa_diag = ones(1,6);
% How kappa = grad sym(psi) is evaluated. All three are computed and
% printed; this one is plotted and used in W:
%   'polyfit' : area-weighted least-squares polynomial fit of sym(psi)
%               of total degree kappaPolyDeg, differentiated analytically.
%               Grid- and mesh-independent; gives the macroscopic curvature,
%               so (L/N)^2*T4 -> 0 like 1/N^2.  RECOMMENDED with mode 5.
%   'element' : mean over translations of the exact within-element strain
%               gradient of every solve. Misses the interface jumps. With
%               full-period averaging the smooth part of the fluctuation
%               gradient does not average out (only smooth part + jumps
%               does), so this T4 grows ~N^2: do not use it with mode 5.
%   'gridFD'  : central differences of sym(psi) on the fixed grid
%               (includes the interface jumps as grid-scale spikes, so
%               T4 depends on nGridX/nGridY)
kappaMethod  = 'polyfit';
kappaPolyDeg = 2;
% Optional [alpha_m alpha_h alpha_c alpha_ell]: if given, W is also
% integrated and compared with ALLSE of the heterogeneous solves.
alphas = [];

% ------------------------------------------------------------
% PLOTTING / OUTPUT
% ------------------------------------------------------------
plotCases         = [];     % cases for the vs-N plots; [] -> all cases
fieldPlotCases    = 1;      % cases for the field figures (strain fields are
                            % stored per translation only for these cases)
plotVsNPerCase    = true;   % one 3-panel vs-N figure per plotted case
sharedCLim        = true;   % Stage-1 and sym(psi) share colour limits per component
overlayInclusions = true;   % dashed zeta=0 inclusions + dots at translated centres
plotOneTranslation = false;  % also plot the raw strain of the middle translation
nContour          = 30;
fieldPlotsN       = [];  % N values for which the field figures are made
closeFieldFigs    = true;   % close field figures once saved (many N x cases)
SAVE_FIGS         = false;
OUTDIR = fullfile(pwd, 'MicroMacroStrain_out');   % results of N go to OUTDIR/N<N>

%% ============ SETUP ============

if ~isfolder(OUTDIR), mkdir(OUTDIR); end

Xc = 0.5*Lx_tot;  Yc = 0.5*Ly_tot;
TOL = 1e-6*max(Lx_tot, Ly_tot);

% ---- homogenized stiffness (PBC on one unit cell, independent of N) ----
mat  = struct('E_m',E_m, 'nu_m',nu_m, 'E_i',E_i, 'nu_i',nu_i, 'thickness',THICKNESS);
geomPBC = struct('Rfrac', Rfrac, 'meshFrac', meshFracHom);
[C_hom, Cs33_hom] = homogenize2D_PBC(mat, geomPBC); %#ok<ASGLU>
fprintf('C_hom (Voigt, [S11;S22;S12] = C_hom*[E11;E22;gamma12]) =\n');
disp(C_hom);

% ---- load cases ----
if strcmpi(loadCaseSource, 'library')
    [rawVecs, caseLabels] = buildCaseLibrary(caseLibrary);
else
    rawVecs = loadRawVectors(loadCaseSource, numCasesFromFile, manualRawVectors, sampleFileName);
    caseLabels = arrayfun(@(k) rawLabel(rawVecs(k,:)), (1:size(rawVecs,1))', 'UniformOutput', false);
end
nCases  = size(rawVecs,1);
caseNames = arrayfun(@(k) sprintf('NC_k%04d', k-1), (1:nCases)', 'UniformOutput', false);
fprintf('Load cases:\n');
for ic = 1:nCases, fprintf('  %s  %s\n', caseNames{ic}, caseLabels{ic}); end
Hc = cell(nCases,1);  Gc = cell(nCases,1);
for ic = 1:nCases
    [Hc{ic}, Gc{ic}] = assembleHG(rawVecs(ic,:), strain0, gScaleFactor);
end
fprintf('Loaded %d load case(s) from source ''%s''.\n', nCases, loadCaseSource);

if isempty(plotCases), plotCases = 1:nCases; end
plotCases = plotCases(plotCases >= 1 & plotCases <= nCases);
fieldPlotCases = fieldPlotCases(fieldPlotCases >= 1 & fieldPlotCases <= nCases);

% ---- fixed evaluation grid ----
xg = linspace(0, Lx_tot, nGridX);
yg = linspace(0, Ly_tot, nGridY);
[XG, YG] = meshgrid(xg, yg);                  % nGridY x nGridX
nP = numel(XG);
Wgrid = trapWeights1D(yg(:)) * trapWeights1D(xg(:)).';
Wgrid = Wgrid(:) / (Lx_tot*Ly_tot);           % area-mean weights (sum = 1)

%% ============ STAGE 1: C_hom solve under KUBC ============

meshSizeHom = meshFracHom * min(Lx_tot, Ly_tot);
[nodesH, elems6H] = buildPlainMesh(Lx_tot, Ly_tot, meshSizeHom);
elemDataH = buildElemData(nodesH, elems6H);
K_hom = assembleGlobalK_LE(elems6H, elemDataH, repmat({C_hom}, size(elems6H,1), 1), ...
                           THICKNESS, 2*size(nodesH,1));
bndH  = boundaryNodes(nodesH, Lx_tot, Ly_tot, TOL);
U_hom = solveKUBC(K_hom, nodesH, bndH, Hc, Gc, Xc, Yc);          % ndof x nCases
ALLSE_hom = 0.5*sum(U_hom .* (K_hom*U_hom), 1).';                 % nCases x 1

Bg_hom  = gridStrainOperator(nodesH, elems6H, xg, yg);            % 3nP x ndof
EpsM_FE = reshape(Bg_hom * U_hom, nP, 3, nCases);                 % [e11 e22 e12]

EpsAff = zeros(nP, 3, nCases);     % sym(H + G.Xc) on the grid
for ic = 1:nCases
    EpsAff(:,:,ic) = affineStrain(Hc{ic}, Gc{ic}, XG(:)-Xc, YG(:)-Yc);
end
fprintf('Stage 1 done: %d nodes, %d T6 elements, %d load case(s).\n', ...
    size(nodesH,1), size(elems6H,1), nCases);

%% ============ LOOP OVER N (Stages 2-4) ============

nN = numel(Nlist);
VN = struct('T1', nan(nN,nCases), 'T2', nan(nN,nCases), 'T3', nan(nN,nCases), ...
            'T4', nan(nN,nCases), 'T5', nan(nN,nCases), 'ell2T4', nan(nN,nCases), ...
            'T4all', nan(nN,nCases,3), 'twoU_het', nan(nN,nCases), ...
            'twoU_het_std', nan(nN,nCases), 'twoU_hom', nan(nN,nCases));

for iN = 1:nN
    N = Nlist(iN);
    OUTDIR_N = fullfile(OUTDIR, sprintf('N%d', N));
    if ~isfolder(OUTDIR_N), mkdir(OUTDIR_N); end
    fprintf('\n=================== N = %d  (%d of %d) ===================\n', N, iN, nN);

    %% ============ STAGE 2: translation-averaged strain = sym(psi) ============

    Lx = Lx_tot/N;  Ly = Ly_tot/N;
    L_cell   = min(Lx, Ly);
    R        = Rfrac    * L_cell;
    meshSize = meshFrac * L_cell;

    runList = buildTranslationList(translationMode, Lx, Ly, R, marginFrac, stepFrac, nPeriod);
    Ntrans  = size(runList,1);
    if translationMode == 5
        % distance (in cell units) of every circle edge to the cell lines,
        % which include the outer boundary: a tiny gap means a sliver element
        sh = unique([runList(:,1)/Lx; runList(:,2)/Ly]);
        ed = [0.5 + sh - Rfrac; 0.5 + sh + Rfrac];
        gapMin = min(abs(ed - round(ed)));
        if gapMin < 0.5*meshFrac
            warning(['Mode 5: an inclusion comes within %.3g cell of the outer ' ...
                     'boundary (element size %.3g cell) -- expect sliver elements; ' ...
                     'change nPeriod if meshing fails.'], gapMin, meshFrac);
        end
    end
    fprintf('Stage 2: N=%d, translationMode=%d -> %d translation(s), %d x %d = %d FE solves.\n', ...
        N, translationMode, Ntrans, nCases, Ntrans, nCases*Ntrans);

    EpsSum   = zeros(nP, 3, nCases);
    GradSum  = zeros(nP, 6, nCases);       % element-wise grad of eps^(k,j)
    ALLSE_tr = zeros(nCases, Ntrans);      % strain energy of every solve
    EpsTrans = zeros(nP, 3, numel(fieldPlotCases), Ntrans);   % only for field-plot cases
    transCenters = cell(Ntrans,1);

    % Loop order: translations outside, load cases inside, so each mesh and
    % stiffness matrix is built once and all K load cases are solved with
    % one multi-RHS solve. The average is taken per load case over the
    % translations only, which gives exactly the note's double loop.
    for it = 1:Ntrans
        zeta1 = runList(it,1);  zeta2 = runList(it,2);

        centers = buildInclusionCenters(N, N, Lx, Ly, zeta1, zeta2, Lx_tot, Ly_tot, R);
        transCenters{it} = centers;
        [nodes, elems6, isInc] = buildRUCMesh(Lx_tot, Ly_tot, centers, R, meshSize);
        elemData = buildElemData(nodes, elems6);

        D_m = planeStrainD(E_m, nu_m);  D_i = planeStrainD(E_i, nu_i);
        Delem = repmat({D_m}, size(elems6,1), 1);
        Delem(isInc) = {D_i};

        K   = assembleGlobalK_LE(elems6, elemData, Delem, THICKNESS, 2*size(nodes,1));
        bnd = boundaryNodes(nodes, Lx_tot, Ly_tot, TOL);
        U   = solveKUBC(K, nodes, bnd, Hc, Gc, Xc, Yc);

        [Bg, Bgrad] = gridStrainOperator(nodes, elems6, xg, yg);
        Eps = reshape(Bg * U, nP, 3, nCases);     % eps^(k,j) on the fixed grid

        EpsSum  = EpsSum + Eps;
        GradSum = GradSum + reshape(Bgrad * U, nP, 6, nCases);
        ALLSE_tr(:,it) = 0.5*sum(U .* (K*U), 1).';
        EpsTrans(:,:,:,it) = Eps(:,:,fieldPlotCases);

        fprintf('  translation %2d/%2d  zeta=(%+.4f,%+.4f)  %6d nodes %6d elems\n', ...
            it, Ntrans, zeta1, zeta2, size(nodes,1), size(elems6,1));
    end
    SymPsi = EpsSum / Ntrans;

    %% ============ STAGE 3: Delta = sym(psi) - eps_M ============

    switch lower(macroRef)
        case 'stage1fe'
            EpsMref = EpsM_FE;   refLabel = '\epsilon_M = Stage-1 FE field';
        case 'symh'
            EpsMref = zeros(nP,3,nCases);
            for ic = 1:nCases
                H = Hc{ic};
                EpsMref(:,:,ic) = repmat([H(1,1), H(2,2), 0.5*(H(1,2)+H(2,1))], nP, 1);
            end
            refLabel = '\epsilon_M = sym(H)';
        case 'affine'
            EpsMref = EpsAff;    refLabel = '\epsilon_M = sym(H + G\cdotX^c)';
        otherwise
            error('Unknown macroRef ''%s''.', macroRef);
    end
    Delta = SymPsi - EpsMref;

    %% ============ STAGE 4: gamma, kappa and the energy terms ============

    Gam  = EpsMref - SymPsi;              % relative strain gamma = eps_M - sym(psi)
    ell  = Lx_tot / N;                    % L/N in the curvature term
    Vol  = Lx_tot * Ly_tot * THICKNESS;
    Cw   = C_hom;                         % stiffness used in the energy terms
    ALLSE_het = mean(ALLSE_tr, 2);

    kapNames = {'element', 'gridFD', 'polyfit'};
    Kap = cell(1,3);
    Kap{1} = GradSum / Ntrans;
    Kap{2} = gridGradient(SymPsi, xg, yg);
    Kap{3} = polyGradient(SymPsi, XG, YG, Wgrid, kappaPolyDeg, Xc, Yc, Lx_tot, Ly_tot);
    iKap = find(strcmpi(kapNames, kappaMethod));
    if isempty(iKap), error('Unknown kappaMethod ''%s''.', kappaMethod); end
    Kappa = Kap{iKap};

    T1 = zeros(nCases,1); T2 = T1; T3 = T1; T5 = T1; T4all = zeros(nCases,3);
    for ic = 1:nCases
        T1(ic) = Vol * (Wgrid.' * quadC(EpsMref(:,:,ic), EpsMref(:,:,ic), Cw));
        T2(ic) = Vol * (Wgrid.' * quadC(SymPsi(:,:,ic),  SymPsi(:,:,ic),  Cw));
        T3(ic) = Vol * (Wgrid.' * quadC(Gam(:,:,ic),     Gam(:,:,ic),     Cw));
        T5(ic) = Vol * (Wgrid.' * quadC(SymPsi(:,:,ic),  Gam(:,:,ic),     Cw));
        for m = 1:3
            T4all(ic,m) = Vol * (Wgrid.' * (Kap{m}(:,:,ic).^2 * Akappa_diag(:)));
        end
    end
    T4 = T4all(:,iKap);
    Wint = [];
    if ~isempty(alphas)
        Wint = alphas(1)/2*T2 + alphas(2)/2*T3 + alphas(3)*T5 + alphas(4)/2*ell^2*T4;
    end

    %% ============ CONSOLE SUMMARY ============

    for ic = 1:nCases
        H = Hc{ic};  G = Gc{ic};
        symH = [H(1,1), H(2,2), 0.5*(H(1,2)+H(2,1))];
        bForce = kubcBodyForce(C_hom, G);
        fprintf('---------------------------------------------------------------\n');
        fprintf(' Case %d [%s]  raw = [%s]\n', ic, caseNames{ic}, sprintf('%+.3g ', rawVecs(ic,:)));
        fprintf('   H  = [%+.4e %+.4e; %+.4e %+.4e]\n', H(1,1), H(1,2), H(2,1), H(2,2));
        fprintf('   G1 = [G111 G112 G122] = [%+.4e %+.4e %+.4e]\n', G(1,1,1), G(1,1,2), G(1,2,2));
        fprintf('   G2 = [G211 G212 G222] = [%+.4e %+.4e %+.4e]\n', G(2,1,1), G(2,1,2), G(2,2,2));
        fprintf('   body force needed by the KUBC polynomial in C_hom: b = [%+.4e %+.4e]\n', bForce);
        fprintf('   component order below: [11 22 12(tensor)]\n');
        fprintf('   sym(H)                         = [%+.4e %+.4e %+.4e]\n', symH);
        % KUBC => <eps>_V = sym(H) exactly for every solve; the grid (trapezoid)
        % mean reproduces it up to the quadrature error at material interfaces
        fprintf('   grid mean Stage-1 eps_M  (~sym(H)) = [%+.4e %+.4e %+.4e]\n', Wgrid.'*EpsM_FE(:,:,ic));
        fprintf('   grid mean sym(psi)       (~sym(H)) = [%+.4e %+.4e %+.4e]\n', Wgrid.'*SymPsi(:,:,ic));
        fprintf('   max|Stage-1 - sym(H+G.X)|      = [%.3e %.3e %.3e]\n', max(abs(EpsM_FE(:,:,ic)-EpsAff(:,:,ic)),[],1));
        fprintf('   max|Delta|  (ref: %-8s)     = [%.3e %.3e %.3e]\n', macroRef, max(abs(Delta(:,:,ic)),[],1));
        fprintf('   RMS  Delta                     = [%.3e %.3e %.3e]\n', sqrt(Wgrid.'*(Delta(:,:,ic).^2)));
        fprintf('   energy terms (domain integrals, C = C_hom, L/N = %.4g):\n', ell);
        fprintf('     T1 = int epsM''C epsM        = %.6e\n', T1(ic));
        fprintf('     T2 = int psi''C psi          = %.6e   (T2/T1 = %.4f)\n', T2(ic), T2(ic)/T1(ic));
        fprintf('     T3 = int gamma''C gamma      = %.6e   (T3/T1 = %.4f)\n', T3(ic), T3(ic)/T1(ic));
        fprintf('     T5 = int psi''C gamma        = %.6e   (T5/T1 = %.4f)\n', T5(ic), T5(ic)/T1(ic));
        fprintf('     T4 = int kappa''A kappa      = %.6e   [%s]   (L/N)^2*T4/T1 = %.4f\n', ...
            T4(ic), kappaMethod, ell^2*T4(ic)/T1(ic));
        fprintf('          T4 by method  element / gridFD / polyfit(deg %d) = %.4e / %.4e / %.4e\n', ...
            kappaPolyDeg, T4all(ic,:));
        fprintf('     check T1-(T2+2T5+T3)        = %.2e (relative)\n', (T1(ic)-T2(ic)-2*T5(ic)-T3(ic))/T1(ic));
        fprintf('     T3+T5 = int epsM''C gamma    = %+.3e   (%.2e of T3; -> 0 for macroRef=stage1FE)\n', ...
            T3(ic)+T5(ic), (T3(ic)+T5(ic))/T3(ic));
        fprintf('     2*ALLSE Stage 1 (C_hom)     = %.6e   (T1 grid-quadrature error %.2e)\n', ...
            2*ALLSE_hom(ic), (T1(ic)-2*ALLSE_hom(ic))/(2*ALLSE_hom(ic)));
        fprintf('     2*ALLSE heterogeneous       = %.6e   (mean over %d translations, std %.2e)\n', ...
            2*ALLSE_het(ic), Ntrans, 2*std(ALLSE_tr(ic,:)));
        if ~isempty(Wint)
            fprintf('     W(alphas) integrated        = %.6e   vs ALLSE_het = %.6e\n', Wint(ic), ALLSE_het(ic));
        end
    end
    fprintf('---------------------------------------------------------------\n');

    % ---- one-line-per-case comparison table ----
    fprintf('\n  ENERGY TERMS PER LOAD CASE  (C = C_hom, A_kappa = diag[%s], kappa: %s, L/N = %.4g)\n', ...
        sprintf('%g ', Akappa_diag), kappaMethod, ell);
    fprintf('  %-9s %12s %12s %12s %12s %12s %12s %12s\n', 'Case', 'T1:eCe', 'T2:psiCpsi', ...
        'T3:gCg', 'T5:psiCg', '(L/N)^2*T4', '2ALLSE_hom', '2ALLSE_het');
    for ic = 1:nCases
        fprintf('  %-9s %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n', caseNames{ic}, ...
            T1(ic), T2(ic), T3(ic), T5(ic), ell^2*T4(ic), 2*ALLSE_hom(ic), 2*ALLSE_het(ic));
    end
    fprintf('\n');

    % ---- CSV ----
    hdr = {'Case','LoadCase','T1_epsCeps','T2_psiCpsi','T3_gamCgam','T5_psiCgam', ...
           'T4_kapAkap','ell2_T4','T4_element','T4_gridFD','T4_polyfit', ...
           'twoALLSE_hom','twoALLSE_het','twoALLSE_het_std','T2_over_T1','T3_over_T1', ...
           'T5_over_T1','ell2T4_over_T1','W_alphas'};
    Wcol = nan(nCases,1);  if ~isempty(Wint), Wcol = Wint; end
    M = [(1:nCases).', T1, T2, T3, T5, T4, ell^2*T4, T4all, 2*ALLSE_hom, 2*ALLSE_het, ...
         2*std(ALLSE_tr,0,2), T2./T1, T3./T1, T5./T1, ell^2*T4./T1, Wcol];
    writeCSV(fullfile(OUTDIR_N, 'EnergyTerms.csv'), hdr, caseNames, M);

    %% ============ PLOTS ============

    cmapSeq = parula(256);
    cmapDiv = divergingMap(256);
    ovl = struct('on', overlayInclusions, 'R', R, ...
        'nominal', buildInclusionCenters(N, N, Lx, Ly, 0, 0, Lx_tot, Ly_tot, R), ...
        'moved', vertcat(transCenters{:}));
    noOvl = ovl; noOvl.on = false;

    if ismember(N, fieldPlotsN), ipList = 1:numel(fieldPlotCases); else, ipList = []; end
    for ip = ipList
        ic  = fieldPlotCases(ip);
        tag = strrep(caseNames{ic}, '_', '\_');
        base = fullfile(OUTDIR_N, caseNames{ic});

        Em = EpsM_FE(:,:,ic);
        Sp = SymPsi(:,:,ic);
        Dl = Delta(:,:,ic);
        zs = max([abs(Em(:)); abs(Sp(:)); realmin]);   % strain scale of this case

        % shared colour limits per tensor component (Stage 1 shows gamma = 2 e12)
        climM = cell(1,3);  climP = cell(1,3);
        for c = 1:3
            if sharedCLim
                lo = min([Em(:,c); Sp(:,c)]);  hi = max([Em(:,c); Sp(:,c)]);
                climP{c} = [lo hi];
            else
                climP{c} = [min(Sp(:,c)) max(Sp(:,c))];
                lo = min(Em(:,c));  hi = max(Em(:,c));
            end
            climM{c} = [lo hi];
        end
        climM{3} = 2*climM{3};

        % ---- Stage 1: macro strain (eps_XX, eps_YY, gamma_XY) ----
        f1 = plotGridTriplet(XG, YG, [Em(:,1), Em(:,2), 2*Em(:,3)], ...
            {'\epsilon_{XX}', '\epsilon_{YY}', '\gamma_{XY}'}, ...
            sprintf('Stage 1: macro strain, C_{hom} + KUBC  --  %s', tag), ...
            ['Stage1_macro_' caseNames{ic}], climM, cmapSeq, nContour, noOvl, zs);

        % ---- Stage 2/3: sym(psi) ----
        f2 = plotGridTriplet(XG, YG, Sp, ...
            {'sym(\psi)_{11}', 'sym(\psi)_{22}', 'sym(\psi)_{12}'}, ...
            sprintf('Stage 2: sym(\\psi), mean over %d translations (N=%d, mode %d)  --  %s', ...
                Ntrans, N, translationMode, tag), ...
            ['Stage2_symPsi_' caseNames{ic}], climP, cmapSeq, nContour, ovl, zs);

        % ---- Stage 3: Delta = sym(psi) - eps_M ----
        climD = cell(1,3);
        for c = 1:3, a = max(abs(Dl(:,c))); climD{c} = [-a a]; end
        f3 = plotGridTriplet(XG, YG, Dl, ...
            {'\Delta_{11}', '\Delta_{22}', '\Delta_{12}'}, ...
            sprintf('Stage 3: \\Delta = sym(\\psi) - \\epsilon_M,  %s  --  %s', refLabel, tag), ...
            ['Stage3_Delta_' caseNames{ic}], climD, cmapDiv, nContour, ovl, zs);

        % ---- midline cuts: every translation, their mean, Stage 1, affine ----
        f4 = plotMidlineCuts(xg, yg, Em, Sp, EpsAff(:,:,ic), squeeze(EpsTrans(:,:,ip,:)), ...
            sprintf('Midline cuts  --  %s', tag), ['Cuts_' caseNames{ic}]);

        % ---- Stage 4: relative strain gamma = eps_M - sym(psi) ----
        Ga = Gam(:,:,ic);
        climG = cell(1,3);
        for c = 1:3, a = max(abs(Ga(:,c))); climG{c} = [-a a]; end
        f6 = plotGridTriplet(XG, YG, Ga, ...
            {'\gamma_{11}', '\gamma_{22}', '\gamma_{12}'}, ...
            sprintf('Stage 4: relative strain \\gamma = \\epsilon_M - sym(\\psi)  --  %s', tag), ...
            ['Stage4_gamma_' caseNames{ic}], climG, cmapDiv, nContour, ovl, zs);

        % ---- Stage 4: curvature kappa = grad sym(psi) ----
        Kc = Kappa(:,:,ic);
        f7 = plotGridTriplet(XG, YG, Kc, ...
            {'\kappa_{xxx} = \partial_X sym(\psi)_{11}', '\kappa_{xxy} = \partial_Y sym(\psi)_{11}', ...
             '\kappa_{yyx} = \partial_X sym(\psi)_{22}', '\kappa_{yyy} = \partial_Y sym(\psi)_{22}', ...
             '\kappa_{xyx} = \partial_X sym(\psi)_{12}', '\kappa_{xyy} = \partial_Y sym(\psi)_{12}'}, ...
            sprintf('Stage 4: \\kappa = \\nabla sym(\\psi)  [%s]  --  %s', kappaMethod, tag), ...
            ['Stage4_kappa_' caseNames{ic}], cell(1,6), cmapSeq, nContour, ovl, ...
            max([abs(Kc(:)); realmin]), [2 3]);

        % ---- Stage 4: energy-term densities and their integrals ----
        dens = [quadC(EpsMref(:,:,ic), EpsMref(:,:,ic), Cw), quadC(Sp, Sp, Cw), ...
                quadC(Ga, Ga, Cw), quadC(Sp, Ga, Cw), ell^2*(Kc.^2*Akappa_diag(:))];
        f8 = plotEnergyTerms(XG, YG, dens, [T1(ic) T2(ic) T3(ic) T5(ic) ell^2*T4(ic)], ...
            2*ALLSE_het(ic), cmapSeq, cmapDiv, nContour, ovl, ...
            sprintf('Stage 4: energy-term densities and integrals  --  %s', tag), ...
            ['Stage4_energy_' caseNames{ic}]);

        figs = [f1 f2 f3 f4 f6 f7 f8];
        names = {'Stage1_macro', 'Stage2_symPsi', 'Stage3_Delta', 'Cuts', ...
                 'Stage4_gamma', 'Stage4_kappa', 'Stage4_energy'};

        % ---- one raw heterogeneous solve (what is being averaged) ----
        if plotOneTranslation
            jm = ceil(Ntrans/2);
            oneOvl = ovl; oneOvl.nominal = transCenters{jm}; oneOvl.moved = zeros(0,2);
            Ej = EpsTrans(:,:,ip,jm);
            f5 = plotGridTriplet(XG, YG, Ej, ...
                {'\epsilon_{11}', '\epsilon_{22}', '\epsilon_{12}'}, ...
                sprintf('Single heterogeneous solve, \\zeta^{(%d)} = (%.3f, %.3f)  --  %s', ...
                    jm, runList(jm,1), runList(jm,2), tag), ...
                ['Single_translation_' caseNames{ic}], {[],[],[]}, cmapSeq, nContour, oneOvl, zs);
            figs(end+1) = f5; names{end+1} = 'SingleTranslation'; %#ok<SAGROW>
        end

        if SAVE_FIGS
            for q = 1:numel(figs)
                print(figs(q), sprintf('%s_%s.png', base, names{q}), '-dpng', '-r150');
            end
            if closeFieldFigs, close(figs); end
        end
    end

    % ---- comparison of the integrated terms across load cases ----
    if nCases > 1 && ismember(N, fieldPlotsN)
        fc = plotEnergyComparison(caseNames, [T1 T2 T3 T5 ell^2*T4 2*ALLSE_het]);
        if SAVE_FIGS
            print(fc, fullfile(OUTDIR_N, 'EnergyTerms_allCases.png'), '-dpng', '-r150');
            if closeFieldFigs, close(fc); end
        end
    end

    %% ============ SAVE RESULTS ============

    gridShape = [nGridY, nGridX];
    res = struct();
    res.xg = xg;  res.yg = yg;  res.Xc = Xc;  res.Yc = Yc;
    res.caseNames = caseNames;  res.rawVecs = rawVecs;  res.H = Hc;  res.G = Gc;
    res.strain0 = strain0;  res.gScaleFactor = gScaleFactor;
    res.C_hom = C_hom;  res.N = N;  res.Rfrac = Rfrac;
    res.translationMode = translationMode;  res.runList = runList;
    res.macroRef = macroRef;
    res.componentOrder = '{11, 22, 12} tensor components (e12 = gamma12/2)';
    res.epsM_stage1 = reshape(EpsM_FE, [gridShape, 3, nCases]);   % ny x nx x 3 x K
    res.epsAffine   = reshape(EpsAff,  [gridShape, 3, nCases]);
    res.symPsi      = reshape(SymPsi,  [gridShape, 3, nCases]);
    res.Delta       = reshape(Delta,   [gridShape, 3, nCases]);
    res.gamma       = reshape(Gam,     [gridShape, 3, nCases]);
    res.kappaOrder  = '[k_xxx k_xxy k_yyx k_yyy k_xyx k_xyy] = d/dX,d/dY of psi11, psi22, psi12';
    res.kappa_element = reshape(Kap{1}, [gridShape, 6, nCases]);
    res.kappa_gridFD  = reshape(Kap{2}, [gridShape, 6, nCases]);
    res.kappa_polyfit = reshape(Kap{3}, [gridShape, 6, nCases]);
    res.kappaMethod = kappaMethod;  res.Akappa_diag = Akappa_diag;  res.ell = ell;
    res.T1 = T1; res.T2 = T2; res.T3 = T3; res.T4 = T4; res.T5 = T5; res.T4all = T4all;
    res.ALLSE_hom = ALLSE_hom;  res.ALLSE_het = ALLSE_het;  res.ALLSE_trans = ALLSE_tr;
    res.alphas = alphas;  res.W_alphas = Wint;
    save(fullfile(OUTDIR_N, 'MicroMacroStrain_results.mat'), '-struct', 'res');
    fprintf('Saved results and figures to %s\n', OUTDIR_N);

    % ---- store the integrated terms of this N ----
    VN.T1(iN,:) = T1.';   VN.T2(iN,:) = T2.';   VN.T3(iN,:) = T3.';
    VN.T4(iN,:) = T4.';   VN.T5(iN,:) = T5.';   VN.ell2T4(iN,:) = (ell^2*T4).';
    VN.T4all(iN,:,:) = reshape(T4all, 1, nCases, 3);
    VN.twoU_het(iN,:) = 2*ALLSE_het.';  VN.twoU_het_std(iN,:) = 2*std(ALLSE_tr,0,2).';
    VN.twoU_hom(iN,:) = 2*ALLSE_hom.';
end

%% ============ BEHAVIOUR OF THE ENERGY TERMS WITH N ============

fprintf('\n  ENERGY TERMS vs N  (C = C_hom, kappa: %s, translationMode %d)\n', kappaMethod, translationMode);
for ic = 1:nCases
    fprintf('  %s   [%s]\n', caseNames{ic}, caseLabels{ic});
    fprintf('    %3s %12s %12s %12s %12s %12s %12s %12s %9s %9s\n', 'N', 'T1', 'T2', 'T3', ...
        'T4', 'T5', '(L/N)^2T4', '2U_het', 'T3/T1', '2U_het/T1');
    for iN = 1:nN
        fprintf('    %3d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %9.4f %9.4f\n', Nlist(iN), ...
            VN.T1(iN,ic), VN.T2(iN,ic), VN.T3(iN,ic), VN.T4(iN,ic), VN.T5(iN,ic), ...
            VN.ell2T4(iN,ic), VN.twoU_het(iN,ic), VN.T3(iN,ic)/VN.T1(iN,ic), ...
            VN.twoU_het(iN,ic)/VN.T1(iN,ic));
    end
    if nN >= 2
        fprintf('    log-log slope vs N:  T3 %+.2f   (T2-T1) %+.2f   (L/N)^2T4 %+.2f   (2U_het-T1) %+.2f\n', ...
            logSlope(Nlist, VN.T3(:,ic)), logSlope(Nlist, VN.T2(:,ic)-VN.T1(:,ic)), ...
            logSlope(Nlist, VN.ell2T4(:,ic)), logSlope(Nlist, VN.twoU_het(:,ic)-VN.T1(:,ic)));
    end
end

% ---- CSV: one row per (case, N) ----
fid = fopen(fullfile(OUTDIR, 'EnergyTerms_vsN.csv'), 'wt');
if fid > 0
    fprintf(fid, 'Case,LoadCase,N,T1,T2,T3,T4,T5,ell2_T4,T4_element,T4_gridFD,T4_polyfit,twoU_hom,twoU_het,twoU_het_std\n');
    for ic = 1:nCases
        for iN = 1:nN
            fprintf(fid, '%d,%s,%d', ic, caseNames{ic}, Nlist(iN));
            fprintf(fid, ',%.10e', VN.T1(iN,ic), VN.T2(iN,ic), VN.T3(iN,ic), VN.T4(iN,ic), ...
                VN.T5(iN,ic), VN.ell2T4(iN,ic), squeeze(VN.T4all(iN,ic,:)), ...
                VN.twoU_hom(iN,ic), VN.twoU_het(iN,ic), VN.twoU_het_std(iN,ic));
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end
VN.Nlist = Nlist;  VN.caseNames = caseNames;  VN.caseLabels = caseLabels;
VN.rawVecs = rawVecs;  VN.kappaMethod = kappaMethod;

% ---- log-log decay rates per case (slope of log(term) vs log(N)) ----
SL = nan(nCases, 5);   % [T3, |T5|, T2-T1, (L/N)^2 T4, 2U_het-T1]
for ic = 1:nCases
    SL(ic,:) = [logSlope(Nlist, VN.T3(:,ic)), logSlope(Nlist, -VN.T5(:,ic)), ...
                logSlope(Nlist, VN.T2(:,ic)-VN.T1(:,ic)), logSlope(Nlist, VN.ell2T4(:,ic)), ...
                logSlope(Nlist, VN.twoU_het(:,ic)-VN.T1(:,ic))];
end
VN.slopes = SL;
VN.slopeOrder = '[T3, -T5, T2-T1, (L/N)^2 T4, 2U_het-T1]';
save(fullfile(OUTDIR, 'EnergyTerms_vsN.mat'), '-struct', 'VN');

iL = nN;   % largest N
fprintf('\n  SUMMARY OVER LOAD CASES  (ratios at N = %d, log-log slopes over N = %s)\n', ...
    Nlist(iL), mat2str(Nlist));
fprintf('  %-9s %-28s %9s %9s %11s %10s | %7s %7s %7s %7s %7s\n', 'Case', 'label', 'T2/T1', ...
    'T3/T1', 'l^2T4/T1', '2U_het/T1', 'sT3', 's|T5|', 'sT2-T1', 'sl^2T4', 'sU-T1');
for ic = 1:nCases
    fprintf('  %-9s %-28s %9.4f %9.5f %11.3e %10.4f | %+7.2f %+7.2f %+7.2f %+7.2f %+7.2f\n', ...
        caseNames{ic}, caseLabels{ic}, VN.T2(iL,ic)/VN.T1(iL,ic), VN.T3(iL,ic)/VN.T1(iL,ic), ...
        VN.ell2T4(iL,ic)/VN.T1(iL,ic), VN.twoU_het(iL,ic)/VN.T1(iL,ic), SL(ic,:));
end
fid = fopen(fullfile(OUTDIR, 'Slopes_vsN.csv'), 'wt');
if fid > 0
    fprintf(fid, 'Case,LoadCase,Label,slope_T3,slope_negT5,slope_T2mT1,slope_ell2T4,slope_2UhetmT1\n');
    for ic = 1:nCases
        fprintf(fid, '%d,%s,"%s"', ic, caseNames{ic}, caseLabels{ic});
        fprintf(fid, ',%.6f', SL(ic,:));  fprintf(fid, '\n');
    end
    fclose(fid);
end

% ---- plots vs N ----
if plotVsNPerCase
    for ip = 1:numel(plotCases)
        ic = plotCases(ip);
        fv = plotTermsVsN(Nlist, VN, ic, sprintf('%s  [%s]', caseNames{ic}, caseLabels{ic}), kappaMethod);
        if SAVE_FIGS
            print(fv, fullfile(OUTDIR, sprintf('%s_EnergyTerms_vsN.png', caseNames{ic})), '-dpng', '-r150');
        end
    end
end
if numel(Nlist) >= 2
    fg = plotRatiosGridAllCases(Nlist, VN, plotCases);
    fo = plotRatiosOverlayAllCases(Nlist, VN, plotCases);
    if SAVE_FIGS
        print(fg, fullfile(OUTDIR, 'AllCases_ratios_vsN_grid.png'),    '-dpng', '-r150');
        print(fo, fullfile(OUTDIR, 'AllCases_ratios_vsN_overlay.png'), '-dpng', '-r150');
    end
end
fprintf('Saved N-sweep tables and plots to %s\n', OUTDIR);


% ============================================================
% ======================= LOCAL FUNCTIONS =====================
% ============================================================

% ---------- load cases (same as RUC_FE_2D_LinearElastic.m) ----------

function [raw, lbl] = buildCaseLibrary(name)
    % raw 9-vectors [H11,H22,H12,G1_11,G1_22,G1_12,G2_11,G2_22,G2_12], unit norm
    iG = struct('G1_11',4, 'G1_22',5, 'G1_12',6, 'G2_11',7, 'G2_22',8, 'G2_12',9);
    sgl = {'G1_11','G1_22','G1_12','G2_11','G2_22','G2_12'};
    raw = zeros(0,9);  lbl = {};
    for k = 1:6
        r = zeros(1,9);  r(iG.(sgl{k})) = 1;
        raw(end+1,:) = r;  lbl{end+1,1} = sgl{k}; %#ok<AGROW>
    end
    switch upper(name)
        case 'G6'
            combos = {};
        case 'G16'
            % {component A, sign A, component B, sign B, label}
            combos = { ...
                'G1_12', +1, 'G2_11', -1, 'bend-X  G1_12-G2_11';       % eps12 = 0, eps11 ~ Y
                'G1_12', +1, 'G2_11', +1, 'shearGrad-X  G1_12+G2_11';  % eps12 ~ X
                'G2_12', +1, 'G1_22', -1, 'bend-Y  G2_12-G1_22';       % eps12 = 0, eps22 ~ X
                'G2_12', +1, 'G1_22', +1, 'shearGrad-Y  G2_12+G1_22';  % eps12 ~ Y
                'G1_11', +1, 'G2_12', +1, 'G1_11+G2_12';               % eps11, eps22 both ~ X
                'G1_11', +1, 'G2_12', -1, 'G1_11-G2_12';
                'G1_12', +1, 'G2_22', +1, 'G1_12+G2_22';               % eps11, eps22 both ~ Y
                'G1_12', +1, 'G2_22', -1, 'G1_12-G2_22';
                'G1_11', +1, 'G2_22', +1, 'G1_11+G2_22'};
            r = zeros(1,9);  r(4:9) = 1/sqrt(6);
            raw(end+1,:) = r;  lbl{end+1,1} = 'all six G (equal)';
        case 'G21'
            combos = {};
            for a = 1:5
                for b = a+1:6
                    combos(end+1,:) = {sgl{a}, +1, sgl{b}, +1, [sgl{a} '+' sgl{b}]}; %#ok<AGROW>
                end
            end
        otherwise
            error('Unknown caseLibrary ''%s'' (use G6, G16 or G21).', name);
    end
    for k = 1:size(combos,1)
        r = zeros(1,9);
        r(iG.(combos{k,1})) = combos{k,2};
        r(iG.(combos{k,3})) = combos{k,4};
        raw(end+1,:) = r / norm(r);  lbl{end+1,1} = combos{k,5}; %#ok<AGROW>
    end
    % keep the 'all six' case last for G16
    if strcmpi(name, 'G16')
        k6 = find(strcmp(lbl, 'all six G (equal)'));
        ord = [setdiff(1:numel(lbl), k6), k6];
        raw = raw(ord,:);  lbl = lbl(ord);
    end
end

function lbl = rawLabel(raw)
    nm = {'H11','H22','H12','G1_11','G1_22','G1_12','G2_11','G2_22','G2_12'};
    nz = find(abs(raw) > 1e-12);
    if isempty(nz), lbl = 'zero'; return; end
    if numel(nz) > 3, lbl = sprintf('%d-component mix', numel(nz)); return; end
    parts = arrayfun(@(k) sprintf('%+.2g %s', raw(k), nm{k}), nz, 'UniformOutput', false);
    lbl = strjoin(parts, ' ');
end

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
        rawVecs = allVecs(1:nCases, :);
    else
        error('Unknown loadCaseSource: %s (use ''file'' or ''manual'').', source);
    end
end

function vecs = readSamplePoints(filePath)
    try
        vecs = readmatrix(filePath, 'CommentStyle', '#');
    catch
        fid = fopen(filePath, 'rt');
        if fid < 0, error('Could not open %s.', filePath); end
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

function [H, G] = assembleHG(raw9, strain0, gScaleFactor)
    vH = raw9(1:3) * strain0;
    vG = raw9(4:9) * (strain0 * gScaleFactor);
    H = zeros(2,2);
    H(1,1) = vH(1); H(2,2) = vH(2); H(1,2) = vH(3); H(2,1) = vH(3);
    G = zeros(2,2,2);   % G(i,j,k), symmetric in (j,k)
    G(1,1,1) = vG(1); G(1,2,2) = vG(2); G(1,1,2) = vG(3); G(1,2,1) = vG(3);
    G(2,1,1) = vG(4); G(2,2,2) = vG(5); G(2,1,2) = vG(6); G(2,2,1) = vG(6);
end

% ---------- kinematics of the KUBC polynomial ----------

function E = affineStrain(H, G, xr, yr)
    % sym(grad u) of u_i = H_ij x_j + 1/2 G_ijk x_j x_k  ->  [e11 e22 e12]
    F11 = H(1,1) + G(1,1,1)*xr + G(1,1,2)*yr;
    F12 = H(1,2) + G(1,2,1)*xr + G(1,2,2)*yr;
    F21 = H(2,1) + G(2,1,1)*xr + G(2,1,2)*yr;
    F22 = H(2,2) + G(2,2,1)*xr + G(2,2,2)*yr;
    E = [F11, F22, 0.5*(F12 + F21)];
end

function b = kubcBodyForce(C, G)
    % body force required for u = H.X + 1/2 G:XX to be in equilibrium in
    % a homogeneous material C (Voigt, engineering shear): b = -div(sigma)
    dEdx = [G(1,1,1); G(2,2,1); G(1,2,1) + G(2,1,1)];   % d[e11 e22 g12]/dx
    dEdy = [G(1,1,2); G(2,2,2); G(1,2,2) + G(2,1,2)];   % d[e11 e22 g12]/dy
    dSdx = C*dEdx;  dSdy = C*dEdy;
    b = -[dSdx(1) + dSdy(3), dSdx(3) + dSdy(2)];
end

% ---------- translation list (same rules as RUC_FE_2D_LinearElastic.m) ----------

function runList = buildTranslationList(mode, Lx, Ly, R, marginFrac, stepFrac, nPeriod)
    z1c = Lx/2 - R - marginFrac*Lx;
    z2c = Ly/2 - R - marginFrac*Ly;
    if (z1c < 0 || z2c < 0) && any(mode == [1 2 3])
        error('Critical translation is negative: inclusion or margin too large.');
    end
    stepX = stepFrac*Lx;  stepY = stepFrac*Ly;
    z1 = -z1c:stepX:z1c;  z2 = -z2c:stepY:z2c;
    switch mode
        case 1
            n = min(numel(z1), numel(z2));
            runList = [z1(1:n).', z2(1:n).'];
        case 2
            [A, B] = ndgrid(z1, z2);
            runList = [A(:), B(:)];
        case 3
            sx = stepX:stepX:z1c;  sy = stepY:stepY:z2c;
            runList = [0, 0; sx.', 0*sx.'; -sx.', 0*sx.'; 0*sy.', sy.'; 0*sy.', -sy.'];
        case 4
            runList = [0, 0];
        case 5
            sp = ((0:nPeriod-1) + 0.5)/nPeriod - 0.5;   % uniform over one period
            [A, B] = ndgrid(sp*Lx, sp*Ly);
            runList = [A(:), B(:)];
        otherwise
            error('Unknown translationMode.');
    end
end

% ---------- FE core (same element, Gauss rule and B-matrix as the RUC code) ----------

function elemData = buildElemData(nodes, elems6)
    [GP_L1, GP_L2, GP_W] = gauss3();
    Nelem = size(elems6,1);
    elemData(Nelem) = struct('dNdX',[],'detJ0',[],'w',[]);
    for e = 1:Nelem
        coords6 = nodes(elems6(e,:),:);
        elemData(e).dNdX  = cell(1,3);
        elemData(e).detJ0 = zeros(1,3);
        elemData(e).w     = GP_W;
        for g = 1:3
            [~, dNdX, detJ] = t6ShapeDeriv(GP_L1(g), GP_L2(g), coords6);
            if detJ <= 0
                error('Element %d has non-positive Jacobian at GP %d.', e, g);
            end
            elemData(e).dNdX{g}  = dNdX;
            elemData(e).detJ0(g) = detJ;
        end
    end
end

function D = planeStrainD(E, nu)
    c = E / ((1+nu)*(1-2*nu));
    D = c * [ 1-nu,  nu,        0; ...
              nu,    1-nu,      0; ...
              0,     0,   (1-2*nu)/2 ];
end

function Ke = elementStiffnessLE(ed, Dmat, THICKNESS)
    Ke = zeros(12,12);
    for g = 1:3
        dNdX = ed.dNdX{g};
        B = zeros(3,12);
        B(1,1:2:end) = dNdX(:,1);
        B(2,2:2:end) = dNdX(:,2);
        B(3,1:2:end) = dNdX(:,2);
        B(3,2:2:end) = dNdX(:,1);
        Ke = Ke + (ed.w(g)*ed.detJ0(g)*THICKNESS) * (B.'*Dmat*B);
    end
end

function K = assembleGlobalK_LE(elems6, elemData, Delem, THICKNESS, ndof)
    Nelem = size(elems6,1);
    Ii = zeros(144*Nelem,1); Jj = Ii; Vv = Ii;
    for e = 1:Nelem
        nE = elems6(e,:);
        gd = zeros(12,1); gd(1:2:end) = 2*nE-1; gd(2:2:end) = 2*nE;
        Ke = elementStiffnessLE(elemData(e), Delem{e}, THICKNESS);
        [a, b] = ndgrid(gd, gd);
        idx = (e-1)*144 + (1:144);
        Ii(idx) = a(:); Jj(idx) = b(:); Vv(idx) = Ke(:);
    end
    K = sparse(Ii, Jj, Vv, ndof, ndof);
end

function bnd = boundaryNodes(nodes, Lx, Ly, TOL)
    x = nodes(:,1); y = nodes(:,2);
    bnd = find(abs(x) < TOL | abs(x-Lx) < TOL | abs(y) < TOL | abs(y-Ly) < TOL);
end

function U = solveKUBC(K, nodes, bnd, Hc, Gc, Xc, Yc)
    % KUBC on every boundary node, all load cases in one multi-RHS solve
    ndof = size(K,1);  nC = numel(Hc);
    xr = nodes(bnd,1) - Xc;  yr = nodes(bnd,2) - Yc;
    U = zeros(ndof, nC);
    for ic = 1:nC
        H = Hc{ic}; G = Gc{ic};
        U(2*bnd-1, ic) = H(1,1)*xr + H(1,2)*yr + 0.5*G(1,1,1)*xr.^2 + G(1,1,2)*xr.*yr + 0.5*G(1,2,2)*yr.^2;
        U(2*bnd,   ic) = H(2,1)*xr + H(2,2)*yr + 0.5*G(2,1,1)*xr.^2 + G(2,1,2)*xr.*yr + 0.5*G(2,2,2)*yr.^2;
    end
    isDir = false(ndof,1);  isDir([2*bnd-1; 2*bnd]) = true;
    fr = find(~isDir);  di = find(isDir);
    U(fr,:) = K(fr,fr) \ (-K(fr,di) * U(di,:));
end

% ---------- strain on the fixed grid ----------

function [Bg, Bgrad] = gridStrainOperator(nodes, elems6, xg, yg)
    % Sparse operator: [e11(all grid pts); e22(...); e12(...)] = Bg * u.
    % Each grid point is located in its T6 element and the strain is the
    % exact element strain there (no nodal smoothing, so the jump across
    % the matrix/inclusion interface is kept).
    % Bgrad (optional): [de11/dX; de11/dY; de22/dX; de22/dY; de12/dX;
    % de12/dY] (6 blocks of nP rows) = Bgrad * u, the within-element strain
    % gradient (constant per straight-sided T6 element).
    nx = numel(xg);  ny = numel(yg);  nP = nx*ny;
    [elemOf, L1of, L2of] = locateGridPoints(nodes, elems6, xg, yg);
    if nargout > 1
        Ge = elementStrainGradients(nodes, elems6);    % 6 x 12 x Nelem
        Gi = zeros(72*nP,1); Gj = Gi; Gv = Gi;
        for p = 1:nP
            nE = elems6(elemOf(p),:);
            dofs = zeros(12,1); dofs(1:2:end) = 2*nE-1; dofs(2:2:end) = 2*nE;
            [rr, cc] = ndgrid((0:5)*nP + p, dofs);
            k = (p-1)*72 + (1:72);
            Gi(k) = rr(:);  Gj(k) = cc(:);
            Gv(k) = reshape(Ge(:,:,elemOf(p)), [], 1);
        end
        Bgrad = sparse(Gi, Gj, Gv, 6*nP, 2*size(nodes,1));
    end
    Ii = zeros(24*nP,1); Jj = Ii; Vv = Ii;
    for p = 1:nP
        nE = elems6(elemOf(p),:);
        [~, dNdX] = t6ShapeDeriv(L1of(p), L2of(p), nodes(nE,:));
        d1 = 2*nE(:)-1;  d2 = 2*nE(:);
        k = (p-1)*24;
        Ii(k+(1:6))   = p;        Jj(k+(1:6))   = d1;        Vv(k+(1:6))   = dNdX(:,1);
        Ii(k+(7:12))  = nP+p;     Jj(k+(7:12))  = d2;        Vv(k+(7:12))  = dNdX(:,2);
        Ii(k+(13:24)) = 2*nP+p;   Jj(k+(13:24)) = [d1; d2];  Vv(k+(13:24)) = 0.5*[dNdX(:,2); dNdX(:,1)];
    end
    Bg = sparse(Ii, Jj, Vv, 3*nP, 2*size(nodes,1));
end

function Ge = elementStrainGradients(nodes, elems6)
    % Strain in a straight-sided T6 is linear in (L1,L2), so its gradient
    % is constant: d eps/dX = [eps1-eps3, eps2-eps3] * inv(A), with eps_m
    % the strain at corner m and A = [x1-x3, x2-x3] the chord Jacobian.
    % Rows: [de11/dX de11/dY de22/dX de22/dY de12/dX de12/dY] (tensor e12).
    Nelem = size(elems6,1);
    Ge = zeros(6,12,Nelem);
    Lc = [1 0; 0 1; 0 0];
    for e = 1:Nelem
        coords6 = nodes(elems6(e,:),:);
        Bm = zeros(3,12,3);
        for m = 1:3
            [~, dNdX] = t6ShapeDeriv(Lc(m,1), Lc(m,2), coords6);
            Bm(1,1:2:end,m) = dNdX(:,1);
            Bm(2,2:2:end,m) = dNdX(:,2);
            Bm(3,1:2:end,m) = 0.5*dNdX(:,2);
            Bm(3,2:2:end,m) = 0.5*dNdX(:,1);
        end
        A = [coords6(1,:)-coords6(3,:); coords6(2,:)-coords6(3,:)].';
        Ai = inv(A);
        D1 = Bm(:,:,1) - Bm(:,:,3);  D2 = Bm(:,:,2) - Bm(:,:,3);
        dX = D1*Ai(1,1) + D2*Ai(2,1);
        dY = D1*Ai(1,2) + D2*Ai(2,2);
        Ge(:,:,e) = [dX(1,:); dY(1,:); dX(2,:); dY(2,:); dX(3,:); dY(3,:)];
    end
end

function K6 = gridGradient(F, xg, yg)
    % central differences of the grid field F (nP x 3 x nC) -> nP x 6 x nC
    nx = numel(xg);  ny = numel(yg);  nC = size(F,3);
    K6 = zeros(size(F,1), 6, nC);
    for ic = 1:nC
        for c = 1:3
            [Zx, Zy] = gradient(reshape(F(:,c,ic), ny, nx), xg(2)-xg(1), yg(2)-yg(1));
            K6(:,2*c-1,ic) = Zx(:);  K6(:,2*c,ic) = Zy(:);
        end
    end
end

function K6 = polyGradient(F, XG, YG, W, deg, Xc, Yc, Lx, Ly)
    % area-weighted LSQ fit of every component with a total-degree-'deg'
    % polynomial in s=(X-Xc)/Lx, t=(Y-Yc)/Ly, differentiated analytically
    s = (XG(:)-Xc)/Lx;  t = (YG(:)-Yc)/Ly;
    [a, b] = meshgrid(0:deg, 0:deg);  keep = a(:)+b(:) <= deg;
    a = a(keep).';  b = b(keep).';
    Phi = s.^a .* t.^b;
    Phx = (a .* s.^max(a-1,0) .* t.^b) / Lx;
    Phy = (b .* s.^a .* t.^max(b-1,0)) / Ly;
    nP = size(F,1);  nC = size(F,3);
    Fm = reshape(F, nP, 3*nC);
    coef = (Phi.' * (W .* Phi)) \ (Phi.' * (W .* Fm));
    Kx = reshape(Phx*coef, nP, 3, nC);  Ky = reshape(Phy*coef, nP, 3, nC);
    K6 = zeros(nP, 6, nC);
    K6(:,1:2:end,:) = Kx;  K6(:,2:2:end,:) = Ky;
end

function w = quadC(a, b, C)
    % pointwise a' C b with tensor strains [e11 e22 e12] -> Voigt [e11 e22 2e12]
    av = [a(:,1), a(:,2), 2*a(:,3)];
    bv = [b(:,1), b(:,2), 2*b(:,3)];
    w = sum((av*C) .* bv, 2);
end

function writeCSV(fname, hdr, names, M)
    fid = fopen(fname, 'wt');
    if fid < 0, warning('Could not write %s.', fname); return; end
    fprintf(fid, '%s', hdr{1});  fprintf(fid, ',%s', hdr{2:end});  fprintf(fid, '\n');
    for r = 1:size(M,1)
        fprintf(fid, '%d,%s', M(r,1), names{r});
        fprintf(fid, ',%.10e', M(r,2:end));
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function [elemOf, L1of, L2of] = locateGridPoints(nodes, elems6, xg, yg)
    % Point location on the corner-node (chord) triangulation, which tiles
    % the rectangle exactly. Points on shared edges go to the element with
    % the largest minimum barycentric coordinate.
    xg = xg(:);  yg = yg(:);
    ny = numel(yg);  nP = numel(xg)*ny;
    elemOf = zeros(nP,1);  L1of = zeros(nP,1);  L2of = zeros(nP,1);
    best = -inf(nP,1);
    tolX = 1e-9 * max(xg(end)-xg(1), yg(end)-yg(1));
    for e = 1:size(elems6,1)
        c = nodes(elems6(e,1:3),:);
        ix = find(xg >= min(c(:,1))-tolX & xg <= max(c(:,1))+tolX);
        iy = find(yg >= min(c(:,2))-tolX & yg <= max(c(:,2))+tolX);
        if isempty(ix) || isempty(iy), continue; end
        [IX, IY] = meshgrid(ix, iy);
        IX = IX(:);  IY = IY(:);
        A = [c(1,1)-c(3,1), c(2,1)-c(3,1); c(1,2)-c(3,2), c(2,2)-c(3,2)];
        L = A \ [xg(IX).' - c(3,1); yg(IY).' - c(3,2)];
        mn = min([L; 1 - L(1,:) - L(2,:)], [], 1).';
        pid = (IX-1)*ny + IY;
        upd = mn > best(pid);
        pid = pid(upd);
        best(pid) = mn(upd);
        elemOf(pid) = e;
        L1of(pid) = L(1,upd).';
        L2of(pid) = L(2,upd).';
    end
    if any(elemOf == 0)
        error('locateGridPoints: %d grid point(s) not inside any element.', nnz(elemOf == 0));
    end
    if min(best) < -1e-6
        warning('locateGridPoints: some grid points lie outside the chord triangulation (min bary = %.2e).', min(best));
    end
end

function w = trapWeights1D(p)
    n = numel(p);  w = zeros(n,1);
    if n == 1, return; end
    w(1) = 0.5*(p(2)-p(1));  w(n) = 0.5*(p(n)-p(n-1));
    for i = 2:n-1, w(i) = 0.5*(p(i+1)-p(i-1)); end
end

% ---------- shape functions / quadrature (identical to the RUC code) ----------

function [N, dNdX, detJ, X, Y] = t6ShapeDeriv(L1,L2,coords6)
    L3 = 1-L1-L2;
    N = [L1*(2*L1-1); L2*(2*L2-1); L3*(2*L3-1); 4*L1*L2; 4*L2*L3; 4*L3*L1];
    dN_dL1 = [4*L1-1; 0; -(4*L3-1); 4*L2; -4*L2; 4*(L3-L1)];
    dN_dL2 = [0; 4*L2-1; -(4*L3-1); 4*L1; 4*(L3-L2); -4*L1];
    x = coords6(:,1); y = coords6(:,2);
    J = [dN_dL1'*x, dN_dL2'*x; dN_dL1'*y, dN_dL2'*y];
    detJ = det(J);
    Jinv = inv(J);
    dNdX = [dN_dL1*Jinv(1,1) + dN_dL2*Jinv(2,1), dN_dL1*Jinv(1,2) + dN_dL2*Jinv(2,2)];
    X = N'*x; Y = N'*y;
end

function [L1w,L2w,W] = gauss3()
    L1w = [2/3, 1/6, 1/6];
    L2w = [1/6, 2/3, 1/6];
    W   = [1/6, 1/6, 1/6];
end

% ---------- geometry / mesh (unchanged from RUC_FE_2D_LinearElastic.m) ----------

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
        error('No inclusion circles intersect the domain -- check N/zeta/geometry.');
    end
end

function [nodes, elems6, isInclusionElem] = buildRUCMesh(Lx_tot,Ly_tot,centers,R,meshSize)
    if exist('createpde','file') ~= 2
        error('buildRUCMesh needs the PDE Toolbox (createpde/geometryFromEdges/generateMesh).');
    end
    Nc = size(centers,1);
    gd = [3;4; 0;Lx_tot;Lx_tot;0; 0;0;Ly_tot;Ly_tot];
    names = {'R1'};  sf = 'R1';
    for k = 1:Nc
        gd = [gd, [1; centers(k,1); centers(k,2); R; 0;0;0;0;0;0]]; %#ok<AGROW>
        names{end+1} = sprintf('C%d',k); %#ok<AGROW>
        sf = [sf '+' names{end}]; %#ok<AGROW>
    end
    sf = ['(' sf ')*R1'];          % clip inclusions crossing the outer boundary
    dl = decsg(gd, sf, char(names)');
    model = createpde(1);
    geometryFromEdges(model, dl);
    generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'quadratic');
    nodes  = model.Mesh.Nodes';
    elems6 = model.Mesh.Elements';
    Nelem = size(elems6,1);
    isInclusionElem = false(Nelem,1);
    for e = 1:Nelem
        cn = elems6(e,1:3);
        cx = mean(nodes(cn,1)); cy = mean(nodes(cn,2));
        isInclusionElem(e) = any((centers(:,1)-cx).^2 + (centers(:,2)-cy).^2 <= R^2);
    end
end

function [nodes, elems6] = buildPlainMesh(Lx_tot, Ly_tot, meshSize)
    if exist('createpde','file') ~= 2
        error('buildPlainMesh needs the PDE Toolbox (createpde/geometryFromEdges/generateMesh).');
    end
    gd = [3;4; 0;Lx_tot;Lx_tot;0; 0;0;Ly_tot;Ly_tot];
    dl = decsg(gd, 'R1', char({'R1'})');
    model = createpde(1);
    geometryFromEdges(model, dl);
    generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'quadratic');
    nodes  = model.Mesh.Nodes';
    elems6 = model.Mesh.Elements';
end

% ---------- plotting ----------

function fh = plotGridTriplet(XG, YG, F, titles, supTitle, figName, clims, cmap, nContour, ovl, zscale, layout)
    % zscale: strain scale of the load case; a component whose range is
    % below 1e-9*zscale is round-off and is drawn as a constant field.
    % layout: [rows cols] of subplots (default [1 size(F,2)]).
    if nargin < 12 || isempty(layout), layout = [1 size(F,2)]; end
    fh = figure('Color','w', 'Name', figName, 'Position', [60 60 1500 60+380*layout(1)]);
    for c = 1:size(F,2)
        ax = subplot(layout(1), layout(2), c);
        Z = reshape(F(:,c), size(XG));
        zr = max(Z(:)) - min(Z(:));
        if zr <= 1e-9*zscale
            % (near-)constant field: contourf cannot draw it
            imagesc(XG(1,:), YG(:,1), Z);  set(ax, 'YDir', 'normal');
            v = mean(Z(:));  d = max(1e-3*abs(v), 1e-9*zscale + eps);
            caxis(ax, [v-d, v+d]);
        else
            contourf(XG, YG, Z, nContour, 'LineColor', 'none');
            if ~isempty(clims{c}) && diff(clims{c}) > 0, caxis(ax, clims{c}); end
        end
        colormap(ax, cmap);  colorbar(ax);
        axis(ax, 'equal');  xlim(ax, [XG(1,1) XG(1,end)]);  ylim(ax, [YG(1,1) YG(end,1)]);
        xlabel(ax, 'X');  ylabel(ax, 'Y');
        title(ax, titles{c});
        if ovl.on
            hold(ax, 'on');
            t = linspace(0, 2*pi, 120);
            for k = 1:size(ovl.nominal,1)
                plot(ax, ovl.nominal(k,1) + ovl.R*cos(t), ovl.nominal(k,2) + ovl.R*sin(t), ...
                    'k--', 'LineWidth', 0.8);
            end
            if ~isempty(ovl.moved)
                plot(ax, ovl.moved(:,1), ovl.moved(:,2), 'k.', 'MarkerSize', 8);
            end
            hold(ax, 'off');
        end
    end
    sgtitle(fh, supTitle);
end

function fh = plotMidlineCuts(xg, yg, Em, Sp, Ea, Etr, supTitle, figName)
    % Etr : nP x 3 x Ntrans (every translation of this load case)
    nx = numel(xg);  ny = numel(yg);
    iyM = ceil(ny/2);  ixM = ceil(nx/2);
    if ismatrix(Etr), Etr = reshape(Etr, size(Etr,1), size(Etr,2), 1); end
    nT = size(Etr,3);
    comp = {'11','22','12'};
    fh = figure('Color','w', 'Name', figName, 'Position', [80 60 1500 780]);
    for c = 1:3
        for row = 1:2
            ax = subplot(2,3,(row-1)*3+c);  hold(ax,'on');  grid(ax,'on');
            if row == 1
                s = xg;  pick = @(v) v(iyM,:);  xl = sprintf('X   (Y = %.3f)', yg(iyM));
            else
                s = yg;  pick = @(v) v(:,ixM).';  xl = sprintf('Y   (X = %.3f)', xg(ixM));
            end
            hT = [];
            for j = 1:nT
                hT = plot(ax, s, pick(reshape(Etr(:,c,j), ny, nx)), '-', ...
                    'Color', [0.72 0.72 0.72], 'LineWidth', 0.8);
            end
            hP = plot(ax, s, pick(reshape(Sp(:,c), ny, nx)), 'k-',  'LineWidth', 2.0);
            hM = plot(ax, s, pick(reshape(Em(:,c), ny, nx)), 'b--', 'LineWidth', 1.6);
            hA = plot(ax, s, pick(reshape(Ea(:,c), ny, nx)), 'r:',  'LineWidth', 1.6);
            xlabel(ax, xl);  ylabel(ax, ['\epsilon_{' comp{c} '}  (tensor)']);
            if row == 1 && c == 1
                legend(ax, [hT hP hM hA], {'\epsilon^{(k,j)} single translations', ...
                    'sym(\psi) = mean', '\epsilon_M Stage 1 (C_{hom})', 'sym(H + G\cdotX^c)'}, ...
                    'Location', 'best');
            end
            hold(ax,'off');
        end
    end
    sgtitle(fh, supTitle);
end

function fh = plotEnergyTerms(XG, YG, dens, Tint, E_het, cmapSeq, cmapDiv, nContour, ovl, supTitle, figName)
    % dens: nP x 5 densities [eCe, psiCpsi, gCg, psiCg, (L/N)^2 kAk]
    ttl = {'\epsilon_M^T C \epsilon_M', 'sym(\psi)^T C sym(\psi)', '\gamma^T C \gamma', ...
           'sym(\psi)^T C \gamma', '(L/N)^2 \kappa^T A_\kappa \kappa'};
    fh = figure('Color','w', 'Name', figName, 'Position', [60 60 1500 820]);
    t = linspace(0, 2*pi, 120);
    for c = 1:5
        ax = subplot(2,3,c);
        Z = reshape(dens(:,c), size(XG));
        if max(Z(:)) - min(Z(:)) > 1e-12*max(abs(dens(:)))
            contourf(XG, YG, Z, nContour, 'LineColor', 'none');
        else
            imagesc(XG(1,:), YG(:,1), Z);  set(ax, 'YDir', 'normal');
        end
        if c == 4
            a = max(abs(Z(:)));  if a > 0, caxis(ax, [-a a]); end
            colormap(ax, cmapDiv);
        else
            colormap(ax, cmapSeq);
        end
        colorbar(ax);  axis(ax, 'equal');
        xlim(ax, [XG(1,1) XG(1,end)]);  ylim(ax, [YG(1,1) YG(end,1)]);
        title(ax, sprintf('%s   (\\int = %.3e)', ttl{c}, Tint(c)));
        if ovl.on
            hold(ax, 'on');
            for k = 1:size(ovl.nominal,1)
                plot(ax, ovl.nominal(k,1) + ovl.R*cos(t), ovl.nominal(k,2) + ovl.R*sin(t), 'k--', 'LineWidth', 0.8);
            end
            hold(ax, 'off');
        end
    end
    ax = subplot(2,3,6);
    vals = [Tint(:); E_het];
    bh = bar(ax, 1:6, vals);  set(bh, 'FaceColor', [0.30 0.45 0.70]);
    set(ax, 'XTick', 1:6, 'XTickLabel', {'T1','T2','T3','T5','\ell^2T4','2U_{het}'});
    ylabel(ax, 'domain integral');  grid(ax, 'on');
    title(ax, {'T1 = T2 + 2T5 + T3,  \ell = L/N', 'U_{het} = ALLSE of the heterogeneous solves'});
    sgtitle(fh, supTitle);
end

function fh = plotEnergyComparison(caseNames, V)
    % V: nCases x 6 = [T1 T2 T3 T5 (L/N)^2T4 2ALLSE_het]
    lbl = {'T1: \epsilon_M^TC\epsilon_M', 'T2: \psi^TC\psi', 'T3: \gamma^TC\gamma', ...
           'T5: \psi^TC\gamma', '(L/N)^2T4: \kappa^TA\kappa', '2ALLSE_{het}'};
    col = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125; ...
           0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933];
    nC = size(V,1);
    xl = strrep(caseNames, '_', '\_');
    fh = figure('Color','w', 'Name', 'EnergyTerms_allCases', 'Position', [60 60 1500 820]);
    ax = subplot(2,1,1);
    bh = bar(ax, 1:nC, V, 'grouped');
    for k = 1:numel(bh), set(bh(k), 'FaceColor', col(k,:)); end
    grid(ax, 'on');  ylabel(ax, 'domain integral');
    legend(ax, lbl, 'Location', 'bestoutside');
    set(ax, 'XTick', 1:nC, 'XTickLabel', xl);
    title(ax, 'energy terms per load case');
    ax = subplot(2,1,2);
    rat = abs(V(:,2:end)) ./ V(:,1);
    lo  = 10^floor(log10(min([rat(rat > 0); 1])));
    bh = bar(ax, 1:nC, rat, 'grouped');
    for k = 1:numel(bh), set(bh(k), 'FaceColor', col(k+1,:), 'BaseValue', lo); end
    set(ax, 'YScale', 'log');  ylim(ax, [lo, 3*max(rat(:))]);
    grid(ax, 'on');  ylabel(ax, '|term| / T1  (log)');
    lbl2 = lbl(2:end);  lbl2{3} = '|T5|: |\psi^TC\gamma|';
    legend(ax, lbl2, 'Location', 'bestoutside');
    set(ax, 'XTick', 1:nC, 'XTickLabel', xl);
    title(ax, 'normalised by T1 = \int \epsilon_M^T C \epsilon_M');
end

function p = logSlope(N, v)
    % least-squares slope of log|v| vs log N (NaN if v changes sign or is 0)
    N = N(:);  v = v(:);
    if any(v <= 0) && any(v >= 0), p = NaN; return; end
    c = [log(N), ones(size(N))] \ log(abs(v));
    p = c(1);
end

function fh = plotTermsVsN(Nl, VN, ic, caseName, kappaMethod)
    Nl = Nl(:);
    T1 = VN.T1(:,ic);  T2 = VN.T2(:,ic);  T3 = VN.T3(:,ic);  T4 = VN.T4(:,ic);
    T5 = VN.T5(:,ic);  L4 = VN.ell2T4(:,ic);  U = VN.twoU_het(:,ic);  Us = VN.twoU_het_std(:,ic);
    col = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125; ...
           0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933; 0.3 0.3 0.3];
    fh = figure('Color','w', 'Name', ['EnergyTerms_vsN_' caseName], 'Position', [40 80 1700 520]);

    % (a) the terms themselves
    ax = subplot(1,3,1);  hold(ax,'on');  grid(ax,'on');
    plot(ax, Nl, T1, '-o',  'Color', col(1,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(1,:));
    plot(ax, Nl, T2, '-s',  'Color', col(2,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(2,:));
    plot(ax, Nl, T3, '-^',  'Color', col(3,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(3,:));
    plot(ax, Nl, T5, '-v',  'Color', col(4,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(4,:));
    plot(ax, Nl, L4, '-d',  'Color', col(5,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(5,:));
    he = errorbar(ax, Nl, U, Us);
    set(he, 'Color', col(6,:), 'LineWidth', 1.6, 'Marker', 'p', 'MarkerFaceColor', col(6,:));
    xlabel(ax, 'N (RUCs per side)');  ylabel(ax, 'domain integral');
    set(ax, 'XTick', Nl);  xlim(ax, [min(Nl)-0.3, max(Nl)+0.3]);
    legend(ax, {'T1 = \int\epsilon_M^TC\epsilon_M', 'T2 = \intsym(\psi)^TCsym(\psi)', ...
        'T3 = \int\gamma^TC\gamma', 'T5 = \intsym(\psi)^TC\gamma', ...
        '(L/N)^2 T4 = (L/N)^2\int\kappa^TA_\kappa\kappa', '2U_{het} (\pm std over \zeta)'}, ...
        'Location', 'best');
    title(ax, '(a) energy terms');

    % (b) quantities that should tend to 1
    ax = subplot(1,3,2);  hold(ax,'on');  grid(ax,'on');
    plot(ax, Nl, T2./T1, '-s', 'Color', col(2,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(2,:));
    plot(ax, Nl, U./T1,  '-p', 'Color', col(6,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(6,:));
    plot(ax, Nl, (T2 + 2*T5 + T3)./T1, '--', 'Color', col(7,:), 'LineWidth', 1.0);
    plot(ax, [min(Nl)-0.3, max(Nl)+0.3], [1 1], 'k:');
    xlabel(ax, 'N (RUCs per side)');  ylabel(ax, 'ratio to T1');
    set(ax, 'XTick', Nl);  xlim(ax, [min(Nl)-0.3, max(Nl)+0.3]);
    legend(ax, {'T2 / T1  (sym(\psi) \rightarrow \epsilon_M)', '2U_{het} / T1', ...
        '(T2+2T5+T3)/T1  (identity = 1)'}, 'Location', 'best');
    title(ax, '(b) should tend to 1');

    % (c) quantities that should vanish, log-log with reference slopes
    ax = subplot(1,3,3);  hold(ax,'on');  grid(ax,'on');
    h = [];  lbl = {};
    h(end+1) = plot(ax, Nl, T3./T1,       '-^', 'Color', col(3,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(3,:));
    lbl{end+1} = sprintf('T3/T1  (slope %+.2f)', logSlope(Nl, T3));
    h(end+1) = plot(ax, Nl, abs(T5)./T1,  '-v', 'Color', col(4,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(4,:));
    lbl{end+1} = sprintf('|T5|/T1  (slope %+.2f)', logSlope(Nl, T5));
    h(end+1) = plot(ax, Nl, L4./T1,       '-d', 'Color', col(5,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(5,:));
    lbl{end+1} = sprintf('(L/N)^2T4/T1  (slope %+.2f)', logSlope(Nl, L4));
    h(end+1) = plot(ax, Nl, T4./T1,       ':d', 'Color', col(5,:), 'LineWidth', 1.2);
    lbl{end+1} = sprintf('T4/T1 [%s]  (slope %+.2f)', kappaMethod, logSlope(Nl, T4));
    dU = U - T1;
    if all(dU > 0)
        h(end+1) = plot(ax, Nl, dU./T1, '-p', 'Color', col(6,:), 'LineWidth', 1.6, 'MarkerFaceColor', col(6,:));
        lbl{end+1} = sprintf('(2U_{het}-T1)/T1  (slope %+.2f)', logSlope(Nl, dU));
    end
    if numel(Nl) >= 2
        a = T3(1)/T1(1);
        h(end+1) = plot(ax, Nl, a*(Nl/Nl(1)).^-1, 'k--', 'LineWidth', 0.8);  lbl{end+1} = '\propto N^{-1}';
        b = L4(1)/T1(1);
        h(end+1) = plot(ax, Nl, b*(Nl/Nl(1)).^-2, 'k-.', 'LineWidth', 0.8);  lbl{end+1} = '\propto N^{-2}';
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', Nl);
    xlim(ax, [min(Nl)*0.9, max(Nl)*1.1]);
    xlabel(ax, 'N (RUCs per side)');  ylabel(ax, 'ratio to T1');
    legend(ax, h, lbl, 'Location', 'best');
    title(ax, '(c) should vanish (log-log)');

    sgtitle(fh, sprintf('Energy terms vs number of RUCs  --  %s', strrep(caseName, '_', '\_')));
end

function [ratios, names, sty] = vsNRatios(VN, ic)
    % quantities that should vanish with N, normalised by T1
    T1 = VN.T1(:,ic);
    ratios = [VN.T3(:,ic)./T1, -VN.T5(:,ic)./T1, VN.ell2T4(:,ic)./T1, ...
              (VN.twoU_het(:,ic) - T1)./T1];
    names  = {'T3/T1', '-T5/T1', '(L/N)^2T4/T1', '(2U_{het}-T1)/T1'};
    sty    = {'-^', '-v', '-d', '-p'};
end

function fh = plotRatiosGridAllCases(Nl, VN, cases)
    % one log-log tile per load case, same axes everywhere
    Nl = Nl(:);  K = numel(cases);
    nc = ceil(sqrt(K));  nr = ceil(K/nc);
    col = [0.929 0.694 0.125; 0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933];
    allv = [];
    for ic = cases(:).', r = vsNRatios(VN, ic); allv = [allv; r(r > 0)]; end %#ok<AGROW>
    if isempty(allv), allv = [1e-6; 1]; end
    yl = [10^floor(log10(min(allv))), 10^ceil(log10(max(allv)))];
    fh = figure('Color','w', 'Name', 'AllCases_ratios_vsN_grid', 'Position', [20 20 1800 1150]);
    for t = 1:K
        ic = cases(t);
        ax = subplot(nr, nc, t);  hold(ax,'on');  grid(ax,'on');
        [r, nm, sty] = vsNRatios(VN, ic);
        h = [];
        for q = 1:4
            v = r(:,q);  v(v <= 0) = NaN;       % log axes: drop non-positive points
            h(end+1) = plot(ax, Nl, v, sty{q}, 'Color', col(q,:), 'LineWidth', 1.3, ...
                'MarkerSize', 4, 'MarkerFaceColor', col(q,:)); %#ok<AGROW>
        end
        a = r(1,1);  if ~(a > 0), a = yl(2)/10; end
        h(end+1) = plot(ax, Nl, a*(Nl/Nl(1)).^-1, 'k--', 'LineWidth', 0.7);
        b = r(1,3);  if ~(b > 0), b = yl(1)*10; end
        h(end+1) = plot(ax, Nl, b*(Nl/Nl(1)).^-2, 'k-.', 'LineWidth', 0.7);
        set(ax, 'XScale','log', 'YScale','log', 'XTick', Nl, 'FontSize', 7);
        xlim(ax, [Nl(1)*0.9, Nl(end)*1.1]);  ylim(ax, yl);
        title(ax, sprintf('%d: %s', ic, strrep(VN.caseLabels{ic}, '_', '\_')), 'FontSize', 8);
        if t > (nr-1)*nc, xlabel(ax, 'N'); end
        if mod(t-1, nc) == 0, ylabel(ax, 'ratio to T1'); end
        if t == 1
            legend(ax, h, [nm, {'\propto N^{-1}', '\propto N^{-2}'}], ...
                'Location', 'southwest', 'FontSize', 6);
        end
    end
    sgtitle(fh, 'Ratios to T1 vs N (log-log), all load cases');
end

function fh = plotRatiosOverlayAllCases(Nl, VN, cases)
    % one panel per quantity, all load cases overlaid
    Nl = Nl(:);  K = numel(cases);
    cm = 0.85*hsv(K);
    mk = {'o','s','^','v','d','p','h','>','<','*','x','+'};
    ls = {'-','--'};
    fh = figure('Color','w', 'Name', 'AllCases_ratios_vsN_overlay', 'Position', [30 30 1700 1050]);
    [~, nm] = vsNRatios(VN, cases(1));
    for q = 1:4
        ax = subplot(2,2,q);  hold(ax,'on');  grid(ax,'on');
        h = [];  lg = {};  v1 = [];
        for t = 1:K
            ic = cases(t);
            r = vsNRatios(VN, ic);  v = r(:,q);  v(v <= 0) = NaN;
            h(end+1) = plot(ax, Nl, v, [ls{1+mod(floor((t-1)/numel(mk)),2)} mk{1+mod(t-1,numel(mk))}], ...
                'Color', cm(t,:), 'LineWidth', 1.2, 'MarkerSize', 5); %#ok<AGROW>
            lg{end+1} = sprintf('%d: %s', ic, strrep(VN.caseLabels{ic}, '_', '\_')); %#ok<AGROW>
            v1(end+1) = v(1); %#ok<AGROW>
        end
        g = exp(mean(log(v1(v1 > 0))));  if isempty(g) || ~(g > 0), g = 1; end
        h(end+1) = plot(ax, Nl, g*(Nl/Nl(1)).^-1, 'k--', 'LineWidth', 1.4);  lg{end+1} = '\propto N^{-1}';
        h(end+1) = plot(ax, Nl, g*(Nl/Nl(1)).^-2, 'k-.', 'LineWidth', 1.4);  lg{end+1} = '\propto N^{-2}';
        set(ax, 'XScale','log', 'YScale','log', 'XTick', Nl);
        xlim(ax, [Nl(1)*0.9, Nl(end)*1.1]);
        xlabel(ax, 'N (RUCs per side)');  ylabel(ax, nm{q});
        title(ax, nm{q});
        if q == 2, legend(ax, h, lg, 'Location', 'eastoutside', 'FontSize', 7); end
    end
    sgtitle(fh, 'Ratios to T1 vs N (log-log), all load cases overlaid');
end

function cm = divergingMap(n)
    % blue - white - red
    h = floor(n/2);
    t = linspace(0,1,h).';
    lo = [t, t, ones(h,1)];
    hi = [ones(n-h,1), flipud(linspace(0,1,n-h).'), flipud(linspace(0,1,n-h).')];
    cm = [lo; hi];
end