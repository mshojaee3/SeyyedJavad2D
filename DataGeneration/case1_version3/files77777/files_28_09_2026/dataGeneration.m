function D = dataGeneration(cfg)
% DATAGENERATION  FE data for the identification of Cbar, Chat and G.
%
%   D = dataGeneration(cfg)
%
% For every load case k and every N in cfg.Nlist (N x N RUCs), from
% classical 2D plane-strain FE solves under KUBC
%   u_i = H_ij x_j + 1/2 G_ijk x_j x_k ,   x = X - X_c :
%
%   Stage 1 : eps_M      strain of the plain domain with C_hom (N-independent)
%   Stage 2 : psi        mean strain over nPeriod^2 inclusion shifts spanning
%                        one whole unit cell (translation average)
%             2U_het     mean FE strain energy u'Ku over the shifts
%             Sfe        mean of int sigma_ij dOmega      [s11 s22 s12]      (exact Gauss)
%             Mfe        mean of int sigma_ij x_k dOmega  [s11 x1, s11 x2, s22 x1,
%                                                          s22 x2, s12 x1, s12 x2]
%   Stage 3 : gamma = eps_M - psi,  kappa = grad psi
%             cfg.kappaMethod = 'local' (default): kappa = grad eps_M - grad gamma,
%             grad eps_M exact (element gradients of the Stage-1 FE field, no
%             window, N-independent), grad gamma by moving least squares (locally
%             linear, Gaussian weights of width sigma = kappaLocalWidth*ell);
%             'localPlain': moving least squares applied to psi itself;
%             'polyfit': global polynomial
%             kappa order: [k111 k112 k221 k222 k121 k122], k_ijk = d psi_ij / d x_k
%   Stage 4 : pure integrals of the field components (tensor components,
%             psi_12 = tensorial shear), numbered as in the paper:
%     A1  psi11^2    A2  psi11 psi22  A3  psi22^2   A4  psi12^2
%     A5  gam11^2    A6  gam11 gam22  A7  gam22^2   A8  gam12^2
%     A9  psi11 gam11  A10 psi11 gam22  A11 psi22 gam11  A12 psi22 gam22  A13 psi12 gam12
%     A14 k111^2  A15 k112^2  A16 k122^2  A17 k222^2  A18 k121^2  A19 k221^2
%     I1..I3   psi11, psi22, psi12        I4..I6   gam11, gam22, gam12
%     I7..I12  k111, k112, k221, k222, k121, k122
%     I13..I18 psi11 x1, psi11 x2, psi22 x1, psi22 x2, psi12 x1, psi12 x2
%     I19..I24 gam11 x1, gam11 x2, gam22 x1, gam22 x2, gam12 x1, gam12 x2
%     AM  macroscopic references for checking: AM1..AM4 = the integrals of
%         A1..A4 with psi -> eps_M, AM14..AM19 = those of A14..A19 with
%         kappa -> grad eps_M (same gradient method); AM5..AM13 unused (NaN)
%     IM  macroscopic references of the I terms: IM1..3 = int eps_M, IM7..12 =
%         int grad eps_M, IM13..18 = int eps_M x_k; IM4..6, IM19..24 unused (NaN)
%
% REGIONS: everything above is evaluated over the full domain (D.A, D.I,
% D.Sfe, D.Mfe, D.twoU_het, D.twoU_hom) and, if cfg.interiorBand > 0, also over
% the interior [a, L-a]^2 with a = interiorBand*ell (D.int.*): FE energy from a
% stiffness matrix of the interior Gauss points only, FE stress integrals from
% the interior Gauss points, A and I with exact grid-cell weights.
% D.VolN, D.LrN (and D.int.VolN, D.int.LrN) are the region volume and rms
% moment arm per N; D.int.valid marks N with a non-empty interior.
%
% Results are also saved to cfg.outDir/DataGeneration.mat.
% Requires homogenize2D_PBC.m and the PDE Toolbox (mesh generation only).

t0 = tic;
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end

g  = cfg.geom;   ma = cfg.mat;   tr = cfg.trans;
Lx_tot = g.Lx_tot;  Ly_tot = g.Ly_tot;
Nlist  = cfg.Nlist(:).';
nN     = numel(Nlist);
Xc = 0.5*Lx_tot;  Yc = 0.5*Ly_tot;
TOL = 1e-6*max(Lx_tot, Ly_tot);
Vol = Lx_tot*Ly_tot*ma.thickness;

% ---------------- load cases ----------------
switch lower(cfg.cases.source)
    case 'library'
        [rawVecs, caseLabels] = buildCaseLibrary(cfg.cases.library);
    otherwise
        rawVecs = loadRawVectors(cfg.cases.source, cfg.cases.numFromFile, ...
                                 cfg.cases.manual, cfg.cases.sampleFile);
        caseLabels = arrayfun(@(k) rawLabel(rawVecs(k,:)), (1:size(rawVecs,1))', 'UniformOutput', false);
end
nC = size(rawVecs,1);
caseNames = arrayfun(@(k) sprintf('NC_k%04d', k-1), (1:nC)', 'UniformOutput', false);
Hc = cell(nC,1);  Gc = cell(nC,1);
for ic = 1:nC
    [Hc{ic}, Gc{ic}] = assembleHG(rawVecs(ic,:), cfg.load.strain0, cfg.load.gScaleFactor);
end
fprintf('[dataGeneration] %d load cases, N = %s\n', nC, mat2str(Nlist));

% ---------------- homogenized stiffness, kappa metric ----------------
mat = struct('E_m',ma.E_m, 'nu_m',ma.nu_m, 'E_i',ma.E_i, 'nu_i',ma.nu_i, 'thickness',ma.thickness);
C_hom = homogenize2D_PBC(mat, struct('Rfrac', g.Rfrac, 'meshFrac', g.meshFracHom));
D_m = planeStrainD(ma.E_m, ma.nu_m);  D_i = planeStrainD(ma.E_i, ma.nu_i);

% evaluation grid: set per N below (at least gridPerCell points per cell)
gpc = 0;  if isfield(cfg.num, 'gridPerCell'), gpc = cfg.num.gridPerCell; end
band = 0;  if isfield(cfg, 'interiorBand'), band = cfg.interiorBand; end

% ---------------- STAGE 1: C_hom, KUBC ----------------
[nodesH, elems6H] = buildPlainMesh(Lx_tot, Ly_tot, g.meshFracHom*min(Lx_tot, Ly_tot));
edH  = buildElemData(nodesH, elems6H);
DelH = repmat({C_hom}, size(elems6H,1), 1);
K_h  = assembleGlobalK_LE(elems6H, edH, DelH, ma.thickness, 2*size(nodesH,1));
U_h  = solveKUBC(K_h, nodesH, boundaryNodes(nodesH, Lx_tot, Ly_tot, TOL), Hc, Gc, Xc, Yc);
twoU_hom = sum(U_h .* (K_h*U_h), 1).';                                 % nC x 1
fprintf('[dataGeneration] Stage 1 done (%d elements)\n', size(elems6H,1));

% ---------------- output arrays ----------------
D.A   = nan(nN, nC, 19);     % A1..A19
D.AM  = nan(nN, nC, 19);     % macroscopic references (Stage-1 field eps_M)
D.IM  = nan(nN, nC, 24);     % macroscopic references of the I terms
D.I   = nan(nN, nC, 24);     % I1..I24
D.Sfe = nan(nN, nC, 3);
D.Mfe = nan(nN, nC, 6);
D.twoU_het = nan(nN, nC);  D.twoU_het_std = nan(nN, nC);  D.T1check = nan(nN, nC);
D.ell2 = nan(nN, 1);  D.Ntrans = nan(nN, 1);  D.nGrid = nan(nN, 1);
D.VolN = Vol*ones(nN,1);  D.LrN = max(Lx_tot, Ly_tot)/sqrt(12)*ones(nN,1);
D.int = struct('A', nan(nN,nC,19), 'AM', nan(nN,nC,19), 'IM', nan(nN,nC,24), 'I', nan(nN,nC,24), 'Sfe', nan(nN,nC,3), 'Mfe', nan(nN,nC,6), ...
               'twoU_het', nan(nN,nC), 'twoU_het_std', nan(nN,nC), 'twoU_hom', nan(nN,nC), ...
               'VolN', nan(nN,1), 'LrN', nan(nN,1), 'valid', false(nN,1), 'band', band);

% ---------------- loop over N ----------------
for iN = 1:nN
    N  = Nlist(iN);
    Lx = Lx_tot/N;  Ly = Ly_tot/N;  Lc = min(Lx, Ly);

    % ---- evaluation grid of this N: nGrid points, at least gridPerCell per cell ----
    nG = cfg.num.nGrid;  if gpc > 0, nG = max(nG, gpc*N + 1); end
    xg = linspace(0, Lx_tot, nG);  yg = linspace(0, Ly_tot, nG);
    [XG, YG] = meshgrid(xg, yg);  nP = numel(XG);
    Wg = trapWeights1D(yg(:)) * trapWeights1D(xg(:)).';
    Wg = Wg(:) / (Lx_tot*Ly_tot);                % sum(Wg) = 1 (for the polynomial projection)
    Wfull = regionWeights(xg, yg, [0 Lx_tot 0 Ly_tot]) * ma.thickness;   % int f = Wfull'*f
    xr = XG(:) - Xc;  yr = YG(:) - Yc;           % moment arms
    [BgH, BgradH] = gridStrainOperator(nodesH, elems6H, xg, yg, true);
    EpsM  = reshape(BgH * U_h, nP, 3, nC);
    GradM = reshape(BgradH * U_h, nP, 6, nC);    % exact grad eps_M (element gradients, Stage 1)
    D.nGrid(iN) = nG;
    R  = g.Rfrac*Lc;
    runList = buildTranslationList(tr.mode, Lx, Ly, R, tr.marginFrac, tr.stepFrac, tr.nPeriod);
    Ntr = size(runList,1);
    fprintf('[dataGeneration] N = %d: %d shifts ...', N, Ntr);  tN = tic;
    aI = band*Lc;                                        % interior band width
    hasInt = band > 0 && 2*aI < min(Lx_tot, Ly_tot) - 1e-12;
    boxI = [aI, Lx_tot-aI, aI, Ly_tot-aI];
    if hasInt
        [~, ~, KrH] = regionOps(nodesH, elems6H, edH, DelH, ma.thickness, Xc, Yc, boxI);
        D.int.twoU_hom(iN,:) = sum(U_h .* (KrH*U_h), 1);
        Wint = regionWeights(xg, yg, boxI) * ma.thickness;
        SfeSumI = zeros(3, nC);  MfeSumI = zeros(6, nC);  E2I = zeros(nC, Ntr);
    end

    EpsSum = zeros(nP,3,nC);  GradSum = zeros(nP,6,nC);
    E2 = zeros(nC, Ntr);  SfeSum = zeros(3, nC);  MfeSum = zeros(6, nC);
    wantGrad = strcmpi(cfg.kappaMethod, 'element');
    for it = 1:Ntr
        centers = buildInclusionCenters(N, N, Lx, Ly, runList(it,1), runList(it,2), Lx_tot, Ly_tot, R);
        [nodes, elems6, isInc] = buildRUCMesh(Lx_tot, Ly_tot, centers, R, g.meshFrac*Lc);
        ed = buildElemData(nodes, elems6);
        Delem = repmat({D_m}, size(elems6,1), 1);  Delem(isInc) = {D_i};
        K = assembleGlobalK_LE(elems6, ed, Delem, ma.thickness, 2*size(nodes,1));
        U = solveKUBC(K, nodes, boundaryNodes(nodes, Lx_tot, Ly_tot, TOL), Hc, Gc, Xc, Yc);
        [Bg, Bgrad] = gridStrainOperator(nodes, elems6, xg, yg, wantGrad);
        EpsSum = EpsSum + reshape(Bg*U, nP, 3, nC);
        if wantGrad, GradSum = GradSum + reshape(Bgrad*U, nP, 6, nC); end
        E2(:,it) = sum(U .* (K*U), 1).';
        [Qs, Qm] = regionOps(nodes, elems6, ed, Delem, ma.thickness, Xc, Yc, [0 Lx_tot 0 Ly_tot]);
        SfeSum = SfeSum + Qs*U;
        MfeSum = MfeSum + Qm*U;
        if hasInt
            [Qs, Qm, Kr] = regionOps(nodes, elems6, ed, Delem, ma.thickness, Xc, Yc, boxI);
            SfeSumI = SfeSumI + Qs*U;  MfeSumI = MfeSumI + Qm*U;
            E2I(:,it) = sum(U .* (Kr*U), 1).';
        end
    end
    Psi = EpsSum / Ntr;
    Gam = EpsM - Psi;
    switch lower(cfg.kappaMethod)
        case 'local',   Kap = GradM - localGradient(Gam, xg, yg, cfg.kappaLocalWidth*Lc);
        case 'localplain', Kap = localGradient(Psi, xg, yg, cfg.kappaLocalWidth*Lc);
        case 'polyfit', Kap = polyGradient(Psi, XG, YG, Wg, cfg.kappaPolyDeg, Xc, Yc, Lx_tot, Ly_tot);
        case 'element', Kap = GradSum / Ntr;
        case 'gridfd',  Kap = gridGradient(Psi, xg, yg);
        otherwise, error('Unknown kappaMethod ''%s''.', cfg.kappaMethod);
    end
    % gradient of the macroscopic field (for the AM references)
    wLoc = 0.1;  if isfield(cfg, 'kappaLocalWidth'), wLoc = cfg.kappaLocalWidth; end
    switch lower(cfg.kappaMethod)
        case 'polyfit',    KapM = polyGradient(EpsM, XG, YG, Wg, cfg.kappaPolyDeg, Xc, Yc, Lx_tot, Ly_tot);
        case 'localplain', KapM = localGradient(EpsM, xg, yg, wLoc*Lc);
        otherwise,         KapM = GradM;          % exact, N-independent
    end

    for ic = 1:nC
        p = Psi(:,:,ic);  q = Gam(:,:,ic);  k = Kap(:,:,ic);    % tensor components
        FA = [p(:,1).^2, p(:,1).*p(:,2), p(:,2).^2, p(:,3).^2, ...          % A1..A4
              q(:,1).^2, q(:,1).*q(:,2), q(:,2).^2, q(:,3).^2, ...          % A5..A8
              p(:,1).*q(:,1), p(:,1).*q(:,2), p(:,2).*q(:,1), p(:,2).*q(:,2), p(:,3).*q(:,3), ... % A9..A13
              k(:,1).^2, k(:,2).^2, k(:,6).^2, k(:,4).^2, k(:,5).^2, k(:,3).^2];  % A14..A19
        FI = [p, q, k, ...                                                   % I1..I12
              p(:,1).*xr, p(:,1).*yr, p(:,2).*xr, p(:,2).*yr, p(:,3).*xr, p(:,3).*yr, ...  % I13..I18
              q(:,1).*xr, q(:,1).*yr, q(:,2).*xr, q(:,2).*yr, q(:,3).*xr, q(:,3).*yr];     % I19..I24
        e  = EpsM(:,:,ic);  kM = KapM(:,:,ic);
        FM = [e(:,1).^2, e(:,1).*e(:,2), e(:,2).^2, e(:,3).^2, nan(nP,9), ...
              kM(:,1).^2, kM(:,2).^2, kM(:,6).^2, kM(:,4).^2, kM(:,5).^2, kM(:,3).^2];
        FM(:,5:13) = 0;
        D.A(iN,ic,:)  = Wfull.' * FA;
        D.AM(iN,ic,:) = Wfull.' * FM;  D.AM(iN,ic,5:13) = NaN;
        FIM = [e, zeros(nP,3), kM, ...
               e(:,1).*xr, e(:,1).*yr, e(:,2).*xr, e(:,2).*yr, e(:,3).*xr, e(:,3).*yr, zeros(nP,6)];
        D.IM(iN,ic,:) = Wfull.' * FIM;  D.IM(iN,ic,[4:6 19:24]) = NaN;
        D.I(iN,ic,:) = Wfull.' * FI;
        D.T1check(iN,ic) = Wfull.' * quadC(EpsM(:,:,ic), EpsM(:,:,ic), C_hom);
        if hasInt
            D.int.A(iN,ic,:)  = Wint.' * FA;
            D.int.AM(iN,ic,:) = Wint.' * FM;  D.int.AM(iN,ic,5:13) = NaN;
            D.int.IM(iN,ic,:) = Wint.' * FIM;  D.int.IM(iN,ic,[4:6 19:24]) = NaN;
            D.int.I(iN,ic,:) = Wint.' * FI;
        end
    end
    if hasInt
        D.int.Sfe(iN,:,:) = reshape((SfeSumI / Ntr).', 1, nC, 3);
        D.int.Mfe(iN,:,:) = reshape((MfeSumI / Ntr).', 1, nC, 6);
        D.int.twoU_het(iN,:)     = mean(E2I, 2).';
        D.int.twoU_het_std(iN,:) = std(E2I, 0, 2).';
        D.int.VolN(iN) = (Lx_tot - 2*aI)*(Ly_tot - 2*aI)*ma.thickness;
        D.int.LrN(iN)  = max(Lx_tot - 2*aI, Ly_tot - 2*aI)/sqrt(12);
        D.int.valid(iN) = true;
    end
    D.Sfe(iN,:,:) = reshape((SfeSum / Ntr).', 1, nC, 3);
    D.Mfe(iN,:,:) = reshape((MfeSum / Ntr).', 1, nC, 6);
    D.twoU_het(iN,:)     = mean(E2, 2).';
    D.twoU_het_std(iN,:) = std(E2, 0, 2).';
    D.ell2(iN)   = (Lx_tot/N)^2;
    D.Ntrans(iN) = Ntr;
    fprintf(' done (%.0f s)\n', toc(tN));
end

% ---------------- metadata, checks, save ----------------
D.Nlist = Nlist;  D.caseNames = caseNames;  D.caseLabels = caseLabels;  D.rawVecs = rawVecs;
D.twoU_hom = repmat(twoU_hom.', nN, 1);
D.C_hom = C_hom;  D.kappaMethod = cfg.kappaMethod;  D.Vol = Vol;
D.cfg = cfg;
q1 = max(abs(D.T1check(:) - D.twoU_hom(:)) ./ D.twoU_hom(:));
fprintf('[dataGeneration] check: max|int eps_M''C eps_M - 2U_hom|/2U_hom = %.2e (grid quadrature)\n', q1);
% eps_M does not depend on N: its curvature integral must not either (gradient-operator check)
kM = sum(D.AM(:,:,14:19), 3);                         % nN x nC
mk = mean(kM, 1);  useK = mk > 1e-6*max(mk);           % cases with a real macro curvature (G loads)
kv = (max(kM(:,useK),[],1) - min(kM(:,useK),[],1)) ./ mk(useK);
fprintf(['[dataGeneration] check: variation of int |grad eps_M|^2 over N = %.1f%% (max over cases);\n' ...
         '                 it should be ~0 -- if large, reduce cfg.kappaLocalWidth (window smooths macro gradients)\n'], ...
         100*max(kv(isfinite(kv))));
save(fullfile(cfg.outDir, 'DataGeneration.mat'), '-struct', 'D');
fprintf('[dataGeneration] saved DataGeneration.mat (%.1f min)\n', toc(t0)/60);
end


% ============================================================
% ================= local functions: new =====================
% ============================================================

function v = voigt(e)
    % tensor [e11 e22 e12] -> Voigt [e11 e22 2e12]
    v = [e(:,1), e(:,2), 2*e(:,3)];
end

function [Qs, Qm, Kr] = regionOps(nodes, elems6, ed, Delem, th, Xc, Yc, box)
    % exact (3-point Gauss) integrals over the Gauss points inside the box
    % [x0 x1 y0 y1], linear (quadratic) in u:
    %   Qs*u    = int [s11 s22 s12]
    %   Qm*u    = int [s11 x1, s11 x2, s22 x1, s22 x2, s12 x1, s12 x2]  (x from the centre)
    %   u'*Kr*u = int sigma : eps  (= 2 x strain energy of the box)
    [L1, L2] = gauss3();
    Nsh = zeros(6,3);
    for gp = 1:3, Nsh(:,gp) = t6ShapeDeriv(L1(gp), L2(gp), [0 0;1 0;0 1;.5 0;.5 .5;0 .5]); end
    wantK = nargout > 2;
    tol = 1e-12 * max(abs(box));
    Nel = size(elems6,1);
    Is = zeros(36*Nel,1);  Js = Is;  Vs = Is;
    Im = zeros(72*Nel,1);  Jm = Im;  Vm = Im;
    if wantK, Ik = zeros(144*Nel,1);  Jk = Ik;  Vk = Ik; end
    for e = 1:Nel
        nE = elems6(e,:);
        dofs = zeros(1,12);  dofs(1:2:end) = 2*nE-1;  dofs(2:2:end) = 2*nE;
        xe = nodes(nE,1);  ye = nodes(nE,2);
        Se = zeros(3,12);  Me = zeros(6,12);  Ke = zeros(12,12);
        for gp = 1:3
            X = Nsh(:,gp).'*xe;  Y = Nsh(:,gp).'*ye;
            if X < box(1)-tol || X > box(2)+tol || Y < box(3)-tol || Y > box(4)+tol, continue; end
            dNdX = ed(e).dNdX{gp};
            B = zeros(3,12);
            B(1,1:2:end) = dNdX(:,1);  B(2,2:2:end) = dNdX(:,2);
            B(3,1:2:end) = dNdX(:,2);  B(3,2:2:end) = dNdX(:,1);
            wdet = ed(e).w(gp) * ed(e).detJ0(gp) * th;
            DB = Delem{e} * B;
            x = X - Xc;  y = Y - Yc;
            Se = Se + wdet*DB;
            Me = Me + wdet*[DB(1,:)*x; DB(1,:)*y; DB(2,:)*x; DB(2,:)*y; DB(3,:)*x; DB(3,:)*y];
            if wantK, Ke = Ke + wdet*(B.'*DB); end
        end
        [r, c] = ndgrid(1:3, dofs);  k = (e-1)*36 + (1:36);
        Is(k) = r(:);  Js(k) = c(:);  Vs(k) = Se(:);
        [r, c] = ndgrid(1:6, dofs);  k = (e-1)*72 + (1:72);
        Im(k) = r(:);  Jm(k) = c(:);  Vm(k) = Me(:);
        if wantK
            [r, c] = ndgrid(dofs, dofs);  k = (e-1)*144 + (1:144);
            Ik(k) = r(:);  Jk(k) = c(:);  Vk(k) = Ke(:);
        end
    end
    nd = 2*size(nodes,1);
    Qs = sparse(Is, Js, Vs, 3, nd);
    Qm = sparse(Im, Jm, Vm, 6, nd);
    if wantK, Kr = sparse(Ik, Jk, Vk, nd, nd); end
end

function W = regionWeights(xg, yg, box)
    % area weights of the grid points for the box [x0 x1 y0 y1]: length of each
    % grid cell [x_i - h/2, x_i + h/2] inside the box (trapezoid rule for the
    % full domain, exact for fractional box edges)
    wx = cellW(xg(:), box(1), box(2));  wy = cellW(yg(:), box(3), box(4));
    W  = wy * wx.';  W = W(:);
end

function w = cellW(x, a, b)
    n = numel(x);  lo = zeros(n,1);  hi = zeros(n,1);
    lo(2:n) = 0.5*(x(1:n-1) + x(2:n));  lo(1) = x(1);
    hi(1:n-1) = lo(2:n);  hi(n) = x(n);
    w = max(0, min(hi, b) - max(lo, a));
end

function K6 = localGradient(F, xg, yg, sig)
    % grad of every column of F (nP x 3 x nC, meshgrid order) by moving least
    % squares: at each grid point a locally linear fit a + b.(x - x0) with
    % Gaussian weights exp(-|x - x0|^2 / (2 sig^2)), truncated at 3 sig;
    % returns b in the kappa order [d/dx, d/dy] per component. Exact for
    % affine fields, unbiased up to the boundary.
    xg = xg(:);  yg = yg(:);  nx = numel(xg);  ny = numel(yg);
    h = min(xg(2)-xg(1), yg(2)-yg(1));
    if sig < 1.5*h
        warning('localGradient: sigma = %.3g is below 1.5 grid spacings; using 1.5h.', sig);
        sig = 1.5*h;
    end
    [Kx0, Kx1, Kx2] = kernelMats(xg, sig);
    [Ky0, Ky1, Ky2] = kernelMats(yg, sig);
    O   = ones(ny, nx);
    m11 = Ky0*O*Kx0.';  m12 = Ky0*O*Kx1.';  m13 = Ky1*O*Kx0.';
    m22 = Ky0*O*Kx2.';  m23 = Ky1*O*Kx1.';  m33 = Ky2*O*Kx0.';
    c11 = m22.*m33 - m23.^2;   c12 = -(m12.*m33 - m13.*m23);  c13 = m12.*m23 - m13.*m22;
    c22 = m11.*m33 - m13.^2;   c23 = -(m11.*m23 - m12.*m13);  c33 = m11.*m22 - m12.^2;
    dt  = m11.*c11 + m12.*c12 + m13.*c13;
    nC  = size(F,3);  K6 = zeros(size(F,1), 6, nC);
    for ic = 1:nC
        for c = 1:3
            Z  = reshape(F(:,c,ic), ny, nx);
            r1 = Ky0*Z*Kx0.';  r2 = Ky0*Z*Kx1.';  r3 = Ky1*Z*Kx0.';
            bx = (c12.*r1 + c22.*r2 + c23.*r3) ./ dt;
            by = (c13.*r1 + c23.*r2 + c33.*r3) ./ dt;
            K6(:,2*c-1,ic) = bx(:);  K6(:,2*c,ic) = by(:);
        end
    end
end

function [K0, K1, K2] = kernelMats(x, sig)
    % K0(i,j) = g(x_j - x_i), K1 = g*d, K2 = g*d^2, Gaussian truncated at 3 sig
    d  = x(:).' - x(:);
    g  = exp(-d.^2/(2*sig^2)) .* (abs(d) <= 3*sig);
    K0 = g;  K1 = g.*d;  K2 = g.*d.^2;
end


% ============================================================
% ======= local functions reused from MicroMacroStrain_2D_LE22.m =======
% ============================================================

function [raw, lbl] = buildCaseLibrary(name)
    % raw 9-vectors [H11,H22,H12,G1_11,G1_22,G1_12,G2_11,G2_22,G2_12], unit norm
    iG = struct('G1_11',4, 'G1_22',5, 'G1_12',6, 'G2_11',7, 'G2_22',8, 'G2_12',9);
    sgl = {'G1_11','G1_22','G1_12','G2_11','G2_22','G2_12'};
    raw = zeros(0,9);  lbl = {};
    for k = 1:6
        r = zeros(1,9);  r(iG.(sgl{k})) = 1;
        raw(end+1,:) = r;  lbl{end+1,1} = sgl{k}; %#ok<AGROW>
    end
    nameG = upper(name);
    withH = strncmp(nameG, 'H3', 2);          % 'H3G16' etc.: three H cases in front
    if withH, nameG = nameG(3:end); end
    switch nameG
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
            error('Unknown caseLibrary ''%s'' (G6, G16, G21, optionally prefixed H3, e.g. H3G16).', name);
    end
    for k = 1:size(combos,1)
        r = zeros(1,9);
        r(iG.(combos{k,1})) = combos{k,2};
        r(iG.(combos{k,3})) = combos{k,4};
        raw(end+1,:) = r / norm(r);  lbl{end+1,1} = combos{k,5}; %#ok<AGROW>
    end
    % keep the 'all six' case last for G16
    if strcmp(nameG, 'G16')
        k6 = find(strcmp(lbl, 'all six G (equal)'));
        ord = [setdiff(1:numel(lbl), k6), k6];
        raw = raw(ord,:);  lbl = lbl(ord);
    end
    if withH                                  % H11, H22, H12 single cases
        rawH = [eye(3), zeros(3,6)];
        raw = [rawH; raw];  lbl = [{'H11'; 'H22'; 'H12'}; lbl];
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

function [Bg, Bgrad, elemOf] = gridStrainOperator(nodes, elems6, xg, yg, wantGrad)
    % Sparse operator: [e11(all grid pts); e22(...); e12(...)] = Bg * u.
    % Each grid point is located in its T6 element and the strain is the
    % exact element strain there (no nodal smoothing, so the jump across
    % the matrix/inclusion interface is kept).
    % elemOf: element containing each grid point.
    % Bgrad (only if wantGrad): [de11/dX; de11/dY; de22/dX; de22/dY; de12/dX;
    % de12/dY] (6 blocks of nP rows) = Bgrad * u, the within-element strain
    % gradient (constant per straight-sided T6 element).
    nx = numel(xg);  ny = numel(yg);  nP = nx*ny;
    [elemOf, L1of, L2of] = locateGridPoints(nodes, elems6, xg, yg);
    Bgrad = [];
    if wantGrad
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
