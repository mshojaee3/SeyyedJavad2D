function [D, computed] = fullscale_solution(P, which, outDir, force)
% =========================================================================
%  FULL-SCALE SOLUTIONS  (everything that needs the heterogeneous microstructure, in one file)
%  D = fullscale_solution(P, 'ident')   data for the identification       D = fullscale_solution(P, 'test')   the test problem
%  [D, computed] = fullscale_solution(P, which, outDir, force)
% =========================================================================
%  Plate Lx x Ly, clamped at x = 0 (u = 0), traction t(y) on x = Lx, plane strain, quadratic triangles:
%        div(sigma) = 0,   sigma = D(x) eps(u)          ->  K U = F           D(x) = D_matrix in the matrix, D_incl in the inclusions
%
%  (a) HOMOGENEOUS plate (D = C_hom):     eps_macro = sym grad u_hom   (does not depend on N)
%  (b) HETEROGENEOUS plate with N x (Lx N/Ly) cells, one circular inclusion per cell, for every N in NList:
%        2U = F'U (twice the strain energy), mean edge displacement on x = Lx, clamp reaction R = K U - F = [Rx; Ry; M about (0,Ly/2)]
%  (c) psi and stress:  bin averages (1/NB of a cell) of the strain and stress of the heterogeneous plate, smoothed by Legendre fits
%        psiMode 'translation': mean over psiPeriod^2 translated positions of the inclusion lattice  (psi = mean strain)
%        psiMode 'centred'    : the plate with the inclusion centred in every cell only (no translation)
%
%  SAVING / REUSE (outDir):  every piece is saved in its own MAT-file  outDir/hom_plate.mat  and  outDir/het_N<N>.mat  together with the input it
%  was made with.  A piece is loaded instead of recomputed when its input is unchanged (so adding N = 5 to NList computes only N = 5).
%  force = true recomputes everything.  computed = true if at least one piece was computed.  Without outDir nothing is saved.
%
%  OUTPUT  D.hom : twoU, ux, uy, react, epsCoef (= eps_macro), sigCoef             coefficient arrays: (degY+1) x (degX+1) x 3 x nLoads
%          D.het(iN) : N, nx, Lcell, twoU, ux, uy, react, psiCoef, sigCoef            (all quantities have one column per load case)
%          and in both, the FIELDS (components 11, 22, 12 with engineering shear; 4th dimension = load case):
%            nbx, nby                 : bin grid (square bins; nbx x nby = plate)
%            epsBin, sigBin           : strain and stress, bin means BEFORE smoothing  (nby x nbx x 3 x nLoads;  'translation': mean over the translations)
%            epsSmooth, sigSmooth     : the same fields AFTER the Legendre smoothing, at the bin centres (epsSmooth = psi of the heterogeneous plate)
%            mesh                     : FE solution of the CENTRED plate (no translation): nodes, elems (T6), U = [u1; v1; u2; v2; ...] (nDof x nLoads),
%                                       Gauss-point data Xg, Wg, eps, sig (single; nElem x 3 points [x 3 components x nLoads]; only if P.out.rawFields)
%            isInc (D.het only)       : inclusion flag of every element of mesh
% =========================================================================
pl = P.(which);  mat = P.material;  fe = P.fe;  Lx = pl.Lx;  Ly = pl.Ly;  trc = pl.trc;  Ns = pl.NList;
if nargin < 3, outDir = ''; end
if nargin < 4 || isempty(force), force = false; end
reuse = P.out.reuse && ~force && ~isempty(outDir);
vers = P.out.version;
keepGauss = ~isfield(P.out, 'rawFields') || P.out.rawFields;                  % store the Gauss-point strain/stress of the FE solution (large) or only mesh + displacements
matKey = struct('Em', mat.Em, 'num', mat.num, 'Ei', mat.Ei, 'nui', mat.nui, 'radiusRatio', mat.radiusRatio, 'thickness', mat.thickness);   % C_hom is NOT an input of the heterogeneous plate
fprintf('  full-scale (%s): plate %.3g x %.3g, N = %s, loads: %s, psi: %s\n', which, Lx, Ly, mat2str(Ns), strjoin(pl.loads, ', '), pl.psiMode);
computed = false;  D = struct();

%% (a) homogeneous plate
keyH = struct('version', vers, 'Chom', mat.Chom, 'thickness', mat.thickness, 'Lx', Lx, 'Ly', Ly, 'trc', trc, 'psiDeg', fe.psiDeg, 'keepGauss', keepGauss);
[D.hom, c] = piece(outDir, 'hom_plate.mat', keyH, reuse, @() solve_hom(P, which, keepGauss));
computed = computed || c;

%% (b), (c) heterogeneous plates
np = 0;
if strcmpi(pl.psiMode, 'translation'), np = pl.psiPeriod; end
for iN = 1:numel(Ns)
    N = Ns(iN);
    keyN = struct('version', vers, 'material', matKey, 'Lx', Lx, 'Ly', Ly, 'trc', trc, 'N', N, 'psiMode', pl.psiMode, 'psiPeriod', np, ...
                  'meshFrac', fe.meshFrac, 'NB', fe.NB, 'psiDeg', fe.psiDeg, 'keepGauss', keepGauss);
    [Dn, c] = piece(outDir, sprintf('het_N%d.mat', N), keyN, reuse, @() solve_het(P, which, N, keepGauss));
    computed = computed || c;
    D.het(iN) = Dn;
end
end

function [out, c] = piece(outDir, name, key, reuse, fcn)
% one saved piece of the full-scale data (see stage_cache)
if isempty(outDir)
    out = fcn();  c = true;
else
    [out, c] = stage_cache(fullfile(outDir, name), key, reuse, fcn);
end
end

function Dh = solve_hom(P, which, keepGauss)
% homogeneous plate: solution, bin means of strain and stress, Legendre fits
pl = P.(which);  mat = P.material;  fe = P.fe;  Lx = pl.Lx;  Ly = pl.Ly;  trc = pl.trc;
nby = 50;  nbx = round(nby*Lx/Ly);  assert(abs(nbx*Ly/nby - Lx) < 1e-9, 'Lx/Ly must give an integer number of 50 bins.');
[nodes, elems] = mesh_structured(nbx, nby, Lx, Ly);
h = fem(nodes, elems, false(size(elems, 1), 1), mat.Chom, mat.Chom, Lx, Ly, trc, mat.thickness);
[Ec, Sc, Er, Sr] = acc_fit(acc_add([], h, Lx, nbx, nby), fe.psiDeg);
Dh = struct('twoU', h.twoU, 'ux', h.ux, 'uy', h.uy, 'react', h.react, 'epsCoef', Ec, 'sigCoef', Sc, ...
            'nbx', nbx, 'nby', nby, 'epsBin', Er, 'sigBin', Sr, 'epsSmooth', legendre_grid(Ec, nbx, nby), 'sigSmooth', legendre_grid(Sc, nbx, nby), ...
            'mesh', pack_mesh(h, keepGauss));
end

function Dn = solve_het(P, which, N, keepGauss)
% heterogeneous plate with N cells across the height
pl = P.(which);  mat = P.material;  fe = P.fe;  Lx = pl.Lx;  Ly = pl.Ly;  trc = pl.trc;
Dm = plane_strain(mat.Em, mat.num);  Di = plane_strain(mat.Ei, mat.nui);
ny = N;  nx = round(Lx*ny/Ly);
assert(abs(nx - Lx*ny/Ly) < 1e-9, 'n_x = Lx*N/Ly = %g is not an integer (N = %d).', Lx*ny/Ly, N);
Lxc = Lx/nx;  Lyc = Ly/ny;  Lc = min(Lxc, Lyc);  R = mat.radiusRatio*Lc;  Hmax = fe.meshFrac*Lc;
fprintf('   N = %d (n_x = %d, cell %.4f)', N, nx, Lc);
[nodes, elems, isInc] = mesh_het(nx, ny, Lx, Ly, Lxc, Lyc, R, Hmax, 0, 0);
h0 = fem(nodes, elems, isInc, Dm, Di, Lx, Ly, trc, mat.thickness);                 % centred inclusion
if strcmpi(pl.psiMode, 'centred'), shifts = [0 0];
else
    np = pl.psiPeriod;  sp = ((0:np-1) + 0.5)/np - 0.5;  [Z1, Z2] = ndgrid(sp*Lxc, sp*Lyc);  shifts = [Z1(:), Z2(:)];
end
acc = [];
for k = 1:size(shifts, 1)
    if strcmpi(pl.psiMode, 'centred'), h = h0;
    else
        [nodes, elems, isInc2] = mesh_het(nx, ny, Lx, Ly, Lxc, Lyc, R, Hmax, shifts(k, 1), shifts(k, 2));
        h = fem(nodes, elems, isInc2, Dm, Di, Lx, Ly, trc, mat.thickness);
    end
    acc = acc_add(acc, h, Lx, nx*fe.NB, ny*fe.NB);  fprintf('.');
end
[Ec, Sc, Er, Sr] = acc_fit(acc, fe.psiDeg);
nbx = nx*fe.NB;  nby = ny*fe.NB;
Dn = struct('N', N, 'nx', nx, 'Lcell', Lc, 'twoU', h0.twoU, 'ux', h0.ux, 'uy', h0.uy, 'react', h0.react, 'psiCoef', Ec, 'sigCoef', Sc, ...
            'nbx', nbx, 'nby', nby, 'epsBin', Er, 'sigBin', Sr, 'epsSmooth', legendre_grid(Ec, nbx, nby), 'sigSmooth', legendre_grid(Sc, nbx, nby), ...
            'mesh', pack_mesh(h0, keepGauss), 'isInc', isInc);
fprintf('\n');
end

function m = pack_mesh(h, keepGauss)
% FE solution kept for the saved fields (Gauss-point data in single precision; omitted when keepGauss is false: mesh and displacements only)
m = struct('nodes', h.nodes, 'elems', h.elems, 'U', h.U);
if keepGauss, m.Xg = single(h.Xg);  m.Wg = single(h.Wg);  m.eps = single(h.eps);  m.sig = single(h.sig); end
end

function Dm = plane_strain(E, nu)
c = E/((1 + nu)*(1 - 2*nu));
Dm = c*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
end

function [nodes, elems, isInc] = mesh_het(nx, ny, Lx, Ly, Lxc, Lyc, R, Hmax, z1, z2)
% heterogeneous mesh (PDE Toolbox): circles of radius R at the cell centres shifted by (z1, z2), clipped by the plate; isInc = inclusion element
centers = zeros(0, 2);
for i = -1:nx
    for j = -1:ny
        cx = (i + 0.5)*Lxc + z1;  cy = (j + 0.5)*Lyc + z2;
        if ~(cx + R < 0 || cx - R > Lx || cy + R < 0 || cy - R > Ly), centers(end+1, :) = [cx, cy]; end %#ok<AGROW>
    end
end
nInc = size(centers, 1);  gd = [3; 4; 0; Lx; Lx; 0; 0; 0; Ly; Ly];  names = {'R1'};  sf = 'R1';
for c = 1:nInc
    gd = [gd, [1; centers(c, 1); centers(c, 2); R; 0; 0; 0; 0; 0; 0]]; %#ok<AGROW>
    names{end+1} = sprintf('C%d', c);  sf = [sf '+' names{end}]; %#ok<AGROW>
end
[dl, ~] = decsg(gd, ['(' sf ')*R1'], char(names)');
model = createpde(1);  geometryFromEdges(model, dl);
msh = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', 'quadratic');
nodes = msh.Nodes';  elems = msh.Elements';
cen = (nodes(elems(:,1), :) + nodes(elems(:,2), :) + nodes(elems(:,3), :))/3;
isInc = false(size(elems, 1), 1);
for k = 1:nInc, isInc = isInc | ((cen(:,1) - centers(k,1)).^2 + (cen(:,2) - centers(k,2)).^2 <= (R + 1e-7)^2); end
end

function o = fem(nodes, elems, isInc, D0, D1, Lx, Ly, trc, THK)
% plane-strain FE solution, quadratic triangles, 3-point rule:  K U = F  with  D = D0 (matrix) or D1 (elements with isInc)
% o.eps, o.sig: strain and stress at the Gauss points (nElem x 3 points x 3 components x nLoads), o.Xg, o.Wg: points and areas,
% o.twoU = F'U, o.ux, o.uy = mean displacement of the edge x = Lx, o.react = [Rx; Ry; M about (0,Ly/2)] of the clamp,
% o.nodes, o.elems, o.U = the mesh and the nodal displacements [u1; v1; u2; v2; ...] (nDof x nLoads)
rg = [2/3 1/6 1/6];  sg = [1/6 2/3 1/6];  wg = 1/6;  TOL = 1e-6;
nE = size(elems, 1);  nd = 2*size(nodes, 1);  nL = size(trc, 3);
xe = reshape(nodes(elems(:), 1), nE, 6);  ye = reshape(nodes(elems(:), 2), nE, 6);
dofs = zeros(nE, 12);  dofs(:, 1:2:end) = 2*elems - 1;  dofs(:, 2:2:end) = 2*elems;
inc = double(isInc(:));
d11 = D0(1,1) + inc*(D1(1,1) - D0(1,1));  d12 = D0(1,2) + inc*(D1(1,2) - D0(1,2));
d22 = D0(2,2) + inc*(D1(2,2) - D0(2,2));  d33 = D0(3,3) + inc*(D1(3,3) - D0(3,3));
Ke = zeros(nE, 144);  B = cell(3, 3);  Xg = zeros(nE, 3, 2);  Wg = zeros(nE, 3);
for g = 1:3
    r = rg(g);  s = sg(g);  t = 1 - r - s;
    N   = [r*(2*r-1); s*(2*s-1); t*(2*t-1); 4*r*s; 4*s*t; 4*t*r];
    dNr = [4*r - 1; 0; -(4*t - 1); 4*s; -4*s; 4*(t - r)];
    dNs = [0; 4*s - 1; -(4*t - 1); 4*r; 4*(t - s); -4*r];
    J11 = xe*dNr;  J12 = ye*dNr;  J21 = xe*dNs;  J22 = ye*dNs;  detJ = J11.*J22 - J12.*J21;
    dNx = ( J22*dNr.' - J12*dNs.')./(detJ*ones(1, 6));
    dNy = (-J21*dNr.' + J11*dNs.')./(detJ*ones(1, 6));
    B1 = zeros(nE, 12);  B1(:, 1:2:end) = dNx;                              % eps11
    B2 = zeros(nE, 12);  B2(:, 2:2:end) = dNy;                              % eps22
    B3 = zeros(nE, 12);  B3(:, 1:2:end) = dNy;  B3(:, 2:2:end) = dNx;       % gamma12
    w = detJ*wg*THK;  DB1 = d11.*B1 + d12.*B2;  DB2 = d12.*B1 + d22.*B2;  DB3 = d33.*B3;
    for i = 1:12
        cols = (i-1)*12 + (1:12);
        Ke(:, cols) = Ke(:, cols) + (w.*B1(:, i)).*DB1 + (w.*B2(:, i)).*DB2 + (w.*B3(:, i)).*DB3;
    end
    B{g,1} = B1;  B{g,2} = B2;  B{g,3} = B3;  Xg(:, g, 1) = xe*N;  Xg(:, g, 2) = ye*N;  Wg(:, g) = w;
end
iK = zeros(nE, 144);  jK = zeros(nE, 144);
for i = 1:12, iK(:, (i-1)*12 + (1:12)) = dofs(:, i)*ones(1, 12);  jK(:, (i-1)*12 + (1:12)) = dofs; end
K = sparse(iK(:), jK(:), Ke(:), nd, nd);
% clamp x = 0, consistent loads of the traction on x = Lx
left = find(abs(nodes(:, 1)) < TOL);  fixed = unique([2*left - 1; 2*left]);
re = zeros(0, 3);  pairs = [1 4 2; 2 5 3; 3 6 1];
for ed = 1:3
    n1 = elems(:, pairs(ed,1));  nm = elems(:, pairs(ed,2));  n2 = elems(:, pairs(ed,3));
    m = abs(nodes(n1,1) - Lx) < TOL & abs(nodes(nm,1) - Lx) < TOL & abs(nodes(n2,1) - Lx) < TOL;
    re = [re; n1(m), nm(m), n2(m)]; %#ok<AGROW>
end
Le = abs(nodes(re(:,3), 2) - nodes(re(:,1), 2));
gp = [-sqrt(3/5), 0, sqrt(3/5)];  gw = [5/9, 8/9, 5/9];  F = zeros(nd, nL);
for q = 1:size(re, 1)
    nn = re(q, :).';  yy = nodes(nn, 2);
    for ig = 1:3
        xi = gp(ig);  Ns = [xi*(xi-1)/2, 1 - xi^2, xi*(xi+1)/2];  dN = [xi - 0.5, -2*xi, xi + 0.5];
        y = Ns*yy;  jac = abs(dN*yy);  s = 2*y/Ly - 1;  Pl = [1, s, (3*s^2 - 1)/2];
        for l = 1:nL
            F(2*nn-1, l) = F(2*nn-1, l) + Ns.'*(Pl*trc(1, :, l).')*jac*gw(ig)*THK;
            F(2*nn,   l) = F(2*nn,   l) + Ns.'*(Pl*trc(2, :, l).')*jac*gw(ig)*THK;
        end
    end
end
free = setdiff(1:nd, fixed);  U = zeros(nd, nL);  U(free, :) = K(free, free) \ F(free, :);
edgeMean = @(v) sum((Le/6).*(v(re(:,1), :) + 4*v(re(:,2), :) + v(re(:,3), :)), 1)/sum(Le);
Rf = K*U - F;  Rx = Rf(2*left - 1, :);  Ry = Rf(2*left, :);
eps = zeros(nE, 3, 3, nL);  sig = eps;
for l = 1:nL
    ue = U(:, l);  ue = ue(dofs);
    for g = 1:3
        e1 = sum(B{g,1}.*ue, 2);  e2 = sum(B{g,2}.*ue, 2);  e3 = sum(B{g,3}.*ue, 2);
        eps(:, g, 1, l) = e1;  eps(:, g, 2, l) = e2;  eps(:, g, 3, l) = e3;
        sig(:, g, 1, l) = d11.*e1 + d12.*e2;  sig(:, g, 2, l) = d12.*e1 + d22.*e2;  sig(:, g, 3, l) = d33.*e3;
    end
end
o = struct('eps', eps, 'sig', sig, 'Xg', Xg, 'Wg', Wg, 'twoU', sum(F.*U, 1), 'ux', edgeMean(U(1:2:end, :)), 'uy', edgeMean(U(2:2:end, :)), ...
           'react', [sum(Rx, 1); sum(Ry, 1); -sum((nodes(left, 2) - Ly/2).*Rx, 1)], 'nodes', nodes, 'elems', elems, 'U', U);
end

function acc = acc_add(acc, h, Lx, nbx, nby)
% add the area-weighted bin sums of strain and stress of one FE solution (bins: squares of size Lx/nbx)
nL = size(h.eps, 4);
if isempty(acc), acc.W = zeros(nby, nbx, nL);  acc.E = zeros(nby, nbx, 3, nL);  acc.S = acc.E; end
d = Lx/nbx;  ix = min(max(floor(h.Xg(:,:,1)/d) + 1, 1), nbx);  iy = min(max(floor(h.Xg(:,:,2)/d) + 1, 1), nby);
idx = [iy(:), ix(:)];  w = h.Wg(:);
for l = 1:nL
    acc.W(:, :, l) = acc.W(:, :, l) + accumarray(idx, w, [nby, nbx]);
    for k = 1:3
        e = h.eps(:, :, k, l);  s = h.sig(:, :, k, l);
        acc.E(:, :, k, l) = acc.E(:, :, k, l) + accumarray(idx, w.*e(:), [nby, nbx]);
        acc.S(:, :, k, l) = acc.S(:, :, k, l) + accumarray(idx, w.*s(:), [nby, nbx]);
    end
end
end

function [Ec, Sc, Er, Sr] = acc_fit(acc, deg)
% bin means (Er, Sr: BEFORE smoothing) -> least-squares Legendre fits (Ec, Sc: coefficients)  f(x,y) = sum_ij C_ij P_i(y) P_j(x)
Wr = repmat(reshape(acc.W, size(acc.W, 1), size(acc.W, 2), 1, size(acc.W, 3)), [1 1 3 1]);
Er = acc.E./Wr;  Sr = acc.S./Wr;
Ec = legendre_fit(Er, deg);  Sc = legendre_fit(Sr, deg);
end

function C = legendre_fit(Z, deg)
sz = size(Z);  nby = sz(1);  nbx = sz(2);  Z3 = reshape(Z, nby, nbx, []);
s = ((1:nbx)' - 0.5)/nbx*2 - 1;  t = ((1:nby)' - 0.5)/nby*2 - 1;
Px = legendre_basis(s, deg(1));  Py = legendre_basis(t, deg(2));
C3 = zeros(deg(2) + 1, deg(1) + 1, size(Z3, 3));
for k = 1:size(Z3, 3), C3(:, :, k) = (Py \ Z3(:, :, k)) / Px.'; end
C = reshape(C3, [deg(2) + 1, deg(1) + 1, sz(3:end)]);
end
