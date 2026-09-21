function [C, Cs33, info] = homogenize2D_PBC(mat, geom)
% HOMOGENIZE2D_PBC  Homogenized 2D plane-strain stiffness Cij via
% periodic boundary conditions (PBC), for a two-phase RVE: isotropic
% matrix with a periodic square array of circular inclusions.
%
%   [C, Cs33] = homogenize2D_PBC(mat, geom)
%
% INPUTS (both optional; omit or pass [] to use the defaults below)
%   mat.E_m, mat.nu_m   : matrix Young's modulus, Poisson's ratio
%   mat.E_i, mat.nu_i   : inclusion Young's modulus, Poisson's ratio
%   mat.thickness       : out-of-plane thickness (cancels exactly out
%                          of Cij; kept only for API completeness)
%   geom.Lx_tot,Ly_tot  : total RVE domain size          (default 1,1)
%   geom.N              : unit cells per side             (default 1)
%   geom.Rfrac          : inclusion radius / cell size     (default 0.30)
%   geom.meshFrac       : element size / cell size         (default 0.05)
%
% OUTPUT
%   C    : 3x3 homogenized plane-strain stiffness (Voigt convention,
%          gamma12 = 2*E12):
%              [Savg11; Savg22; Savg12] = C * [E11; E22; gamma12]
%   Cs33 : 1x3 row giving the homogenized out-of-plane coupling, from
%          the SAME three unit-strain solves (no extra cost):
%              Savg33 = Cs33 * [E11; E22; gamma12]
%   info : (optional, only computed if requested) struct with the mesh
%          and the three solved displacement fields, for inspection.
%
% METHOD
%   Displacement is split as u(X) = Ebar*(X-Xc) + u*(X) with u*
%   periodic. This gives linear node-pair constraints on opposite
%   edges/corners, u(x+L) - u(x) = Ebar*L, implemented here by direct
%   master/slave elimination (interior nodes + left/bottom edges are
%   master DOFs; right/top edges and three corners are slaved to them;
%   the fourth corner, at the origin, is pinned to remove rigid-body
%   translation). Three canonical unit macroscopic strains (uniaxial-X,
%   uniaxial-Y, shear) are solved on the SAME stiffness matrix; the
%   volume-averaged stress from each solve is one column of C.
%
%   Because the boundary condition is exactly periodic and the
%   microstructure itself is periodic (single circular inclusion per
%   square unit cell), ONE unit cell (N=1) already represents the
%   infinite periodic medium exactly -- no translation ensemble and no
%   RVE-size sweep are needed (unlike KUBC, where both exist only to
%   average out KUBC's own boundary-layer bias).
%
% REQUIRES: MATLAB PDE Toolbox (createpde/geometryFromEdges/
% generateMesh), used only for mesh generation.
%
% NOTE: written and reviewed carefully but not executed here (no
% MATLAB/PDE Toolbox available in this environment). The one place
% this could break is if the automatic mesher does not place matching
% nodes on opposite edges of the unit cell -- pairEdges() below checks
% this explicitly and throws a clear 'MeshNotPeriodic' error rather
% than silently returning wrong numbers, so if you hit that error,
% try a different meshFrac (or report it back and I'll adjust the
% mesher). Please run it and let me know what happens.
%
% EXAMPLE
%   mat.E_m = 70000; mat.nu_m = 0.33;
%   mat.E_i = 3500;  mat.nu_i = 0.33;
%   geom.Rfrac = 0.30; geom.meshFrac = 0.05;
%   [C, Cs33] = homogenize2D_PBC(mat, geom)

if nargin < 1 || isempty(mat),  mat  = struct(); end
if nargin < 2 || isempty(geom), geom = struct(); end

if ~isfield(mat,'E_m'),       mat.E_m = 70000;      end
if ~isfield(mat,'nu_m'),      mat.nu_m = 0.33;       end
if ~isfield(mat,'E_i'),       mat.E_i = 3500;        end
if ~isfield(mat,'nu_i'),      mat.nu_i = 0.33;       end
if ~isfield(mat,'thickness'), mat.thickness = 1.0;   end

if ~isfield(geom,'Lx_tot'),   geom.Lx_tot = 1.0;     end
if ~isfield(geom,'Ly_tot'),   geom.Ly_tot = 1.0;     end
if ~isfield(geom,'N'),        geom.N = 1;            end
if ~isfield(geom,'Rfrac'),    geom.Rfrac = 0.30;     end
if ~isfield(geom,'meshFrac'), geom.meshFrac = 0.05;  end

Lx = geom.Lx_tot; Ly = geom.Ly_tot; N = 1;
THICKNESS = mat.thickness;
L_cell = min(Lx,Ly) / N;
R        = geom.Rfrac    * L_cell;
meshSize = geom.meshFrac * L_cell;
TOL = 1e-6 * max(Lx,Ly);

%% ---- mesh (single periodic unit-cell pattern, tiled N x N) ----
centers = buildPeriodicInclusionCenters(N, Lx, Ly, R);
[nodes, elems6, isInclusionElem] = buildRUCMesh(Lx, Ly, centers, R, meshSize);
Nnodes = size(nodes,1); Nelem = size(elems6,1);
ndof = 2*Nnodes;
V0 = Lx*Ly*THICKNESS;

%% ---- element geometric data (shape derivatives, once per mesh) ----
[GP_L1, GP_L2, GP_W] = gauss3();
elemData(Nelem) = struct('dNdX',[],'detJ0',[],'w',[]);
for e = 1:Nelem
    coords6 = nodes(elems6(e,:),:);
    elemData(e).dNdX = cell(1,3);
    elemData(e).detJ0 = zeros(1,3);
    elemData(e).w = GP_W;
    for g = 1:3
        [~,dNdX,detJ] = t6ShapeDeriv(GP_L1(g),GP_L2(g),coords6);
        if detJ <= 0
            error('homogenize2D_PBC:BadJacobian', ...
                'Element %d has non-positive Jacobian at GP %d.', e, g);
        end
        elemData(e).dNdX{g} = dNdX;
        elemData(e).detJ0(g) = detJ;
    end
end

%% ---- material per element, global stiffness (assembled once) ----
D_m = planeStrainD(mat.E_m, mat.nu_m);
D_i = planeStrainD(mat.E_i, mat.nu_i);
Delem = cell(Nelem,1); nu_e = zeros(Nelem,1);
for e = 1:Nelem
    if isInclusionElem(e), Delem{e} = D_i; nu_e(e) = mat.nu_i;
    else,                  Delem{e} = D_m; nu_e(e) = mat.nu_m;
    end
end
K = assembleGlobalK_LE(elems6, elemData, Delem, THICKNESS, ndof);

%% ---- periodic boundary bookkeeping (geometry-only, once) ----
[leftE,rightE,botE,topE,c00,c10,c01,c11] = classifyPeriodicBoundary(nodes, Lx, Ly, TOL);
[leftE,rightE] = pairEdges(nodes, leftE, rightE, 2, TOL);   % match by Y
[botE,topE]    = pairEdges(nodes, botE,  topE,   1, TOL);   % match by X

%% ---- solve the 3 canonical unit macroscopic-strain cases ----
strainCases = eye(3);   % rows: [E11,E22,gamma12] = [1,0,0]/[0,1,0]/[0,0,1]
C = zeros(3,3);
Cs33 = zeros(1,3);
wantInfo = nargout > 2;
if wantInfo
    info.u = cell(1,3);
    info.strainRecoveryResid = zeros(1,3);
end

for k = 1:3
    e11 = strainCases(k,1); e22 = strainCases(k,2); g12 = strainCases(k,3);
    Ebar = [e11, g12/2; g12/2, e22];

    [T, up] = buildPeriodicReduction(ndof, leftE,rightE,botE,topE, ...
                                      c00,c10,c01,c11, Lx,Ly, Ebar);
    Ku = K*up;
    Kr = T.' * K * T;
    Fr = -T.' * Ku;
    umaster = Kr \ Fr;
    u = T*umaster + up;

    [Savg, Eavg] = domainAverageStressStrain(u, elems6, elemData, Delem, nu_e, THICKNESS, V0);
    C(:,k) = [Savg(1); Savg(2); Savg(4)];
    Cs33(k) = Savg(3);

    resid = norm([Eavg(1);Eavg(2);2*Eavg(4)] - strainCases(k,:).');
    if resid > 1e-8
        warning('homogenize2D_PBC:StrainRecovery', ...
            'Case %d: recovered average strain differs from the prescribed one by %.2e -- check the periodic BC setup.', k, resid);
    end
    if wantInfo
        info.u{k} = u;
        info.strainRecoveryResid(k) = resid;
    end
end

if max(max(abs(C - C.'))) > 1e-6 * max(1, norm(C))
    warning('homogenize2D_PBC:AsymmetricC', ...
        'Homogenized C is not symmetric to within tolerance -- check mesh/BC.');
end

if wantInfo
    info.nodes = nodes; info.elems6 = elems6; info.isInclusionElem = isInclusionElem;
end

end


%% ================= helper functions =================

function centers = buildPeriodicInclusionCenters(N, Lx, Ly, R)
    cellLx = Lx/N; cellLy = Ly/N;
    centers = zeros(0,2);
    for i = -1:N
        for j = -1:N
            cx = (i+0.5)*cellLx;
            cy = (j+0.5)*cellLy;
            if circleIntersectsDomain(cx,cy,R,0,Lx,0,Ly)
                centers(end+1,:) = [cx,cy]; %#ok<AGROW>
            end
        end
    end
end

function tf = circleIntersectsDomain(cx,cy,Rc,xmin,xmax,ymin,ymax)
    tf = ~(cx+Rc < xmin || cx-Rc > xmax || cy+Rc < ymin || cy-Rc > ymax);
end

function [nodes, elems6, isInclusionElem] = buildRUCMesh(Lx_tot,Ly_tot,centers,R,meshSize)
    if exist('createpde','file') ~= 2
        error('homogenize2D_PBC:NoPDEToolbox', ['Needs the Partial Differential ' ...
            'Equation Toolbox (createpde/geometryFromEdges/generateMesh). If ' ...
            'unavailable, supply your own T6 mesh with matching periodic ' ...
            'boundary nodes on opposite edges.']);
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

function [Savg, Eavg] = domainAverageStressStrain(u, elems6, elemData, Delem, nu_e, THICKNESS, V0)
    Nelem = size(elems6,1);
    Sacc = zeros(4,1); Eacc = zeros(4,1);   % [11,22,33,12]
    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        Dmat = Delem{e}; nu = nu_e(e);
        for g = 1:3
            dNdX = elemData(e).dNdX{g};
            B = zeros(3,12);
            for a = 1:6
                B(1, 2*a-1) = dNdX(a,1);
                B(2, 2*a)   = dNdX(a,2);
                B(3, 2*a-1) = dNdX(a,2);
                B(3, 2*a)   = dNdX(a,1);
            end
            eeng = B*ue;
            s = Dmat*eeng;
            s33 = nu*(s(1)+s(2));
            e12 = 0.5*eeng(3);
            wdv = elemData(e).w(g) * elemData(e).detJ0(g) * THICKNESS;
            Sacc = Sacc + wdv*[s(1);s(2);s33;s(3)];
            Eacc = Eacc + wdv*[eeng(1);eeng(2);0;e12];
        end
    end
    Savg = Sacc/V0; Eavg = Eacc/V0;
end

function [leftE,rightE,botE,topE,c00,c10,c01,c11] = classifyPeriodicBoundary(nodes, Lx, Ly, TOL)
    x = nodes(:,1); y = nodes(:,2);
    isLeft  = abs(x-0)  < TOL;
    isRight = abs(x-Lx) < TOL;
    isBot   = abs(y-0)  < TOL;
    isTop   = abs(y-Ly) < TOL;

    c00 = find(isLeft  & isBot);
    c10 = find(isRight & isBot);
    c01 = find(isLeft  & isTop);
    c11 = find(isRight & isTop);
    if ~(isscalar(c00) && isscalar(c10) && isscalar(c01) && isscalar(c11))
        error('homogenize2D_PBC:CornerNotFound', ...
            'Expected exactly one mesh node at each of the four unit-cell corners.');
    end

    leftE  = find(isLeft  & ~isBot & ~isTop);
    rightE = find(isRight & ~isBot & ~isTop);
    botE   = find(isBot   & ~isLeft & ~isRight);
    topE   = find(isTop   & ~isLeft & ~isRight);
end

function [edgeA, edgeB] = pairEdges(nodes, edgeA, edgeB, coordCol, tol)
    cA = nodes(edgeA, coordCol); [cA, ordA] = sort(cA); edgeA = edgeA(ordA);
    cB = nodes(edgeB, coordCol); [cB, ordB] = sort(cB); edgeB = edgeB(ordB);
    if numel(edgeA) ~= numel(edgeB)
        error('homogenize2D_PBC:MeshNotPeriodic', ...
            ['Boundary node counts differ (%d vs %d) on a pair of opposite ' ...
             'edges -- the mesh is not periodicity-conforming. Try a ' ...
             'different meshFrac.'], numel(edgeA), numel(edgeB));
    end
    if max(abs(cA-cB)) > tol
        error('homogenize2D_PBC:MeshNotPeriodic', ...
            ['Paired boundary coordinates differ by up to %.3g (tol=%.3g) ' ...
             'on a pair of opposite edges -- the mesh is not periodicity-' ...
             'conforming.'], max(abs(cA-cB)), tol);
    end
    edgeA = edgeA(:); edgeB = edgeB(:);
end

function [T, up] = buildPeriodicReduction(ndof, leftE,rightE,botE,topE, ...
                                           c00,c10,c01,c11, Lx,Ly, Ebar)
    isFixed = false(ndof,1);
    isFixed(2*c00-1) = true; isFixed(2*c00) = true;

    isSlave = false(ndof,1);
    slaveMasterDof = zeros(ndof,1);   % 0 => slaved directly to a prescribed value
    up = zeros(ndof,1);

    offX  = Ebar*[Lx;0];   % right - left
    offY  = Ebar*[0;Ly];   % top - bottom
    offXY = Ebar*[Lx;Ly];  % corner11 - corner00

    for k = 1:numel(leftE)
        l = leftE(k); r = rightE(k);
        isSlave(2*r-1)=true; slaveMasterDof(2*r-1)=2*l-1; up(2*r-1)=offX(1);
        isSlave(2*r)  =true; slaveMasterDof(2*r)  =2*l;   up(2*r)  =offX(2);
    end
    for k = 1:numel(botE)
        b = botE(k); t = topE(k);
        isSlave(2*t-1)=true; slaveMasterDof(2*t-1)=2*b-1; up(2*t-1)=offY(1);
        isSlave(2*t)  =true; slaveMasterDof(2*t)  =2*b;   up(2*t)  =offY(2);
    end
    isSlave(2*c10-1)=true; up(2*c10-1)=offX(1);
    isSlave(2*c10)  =true; up(2*c10)  =offX(2);
    isSlave(2*c01-1)=true; up(2*c01-1)=offY(1);
    isSlave(2*c01)  =true; up(2*c01)  =offY(2);
    isSlave(2*c11-1)=true; up(2*c11-1)=offXY(1);
    isSlave(2*c11)  =true; up(2*c11)  =offXY(2);

    isMaster = ~isFixed & ~isSlave;
    masterList = find(isMaster);
    nMaster = numel(masterList);
    masterCol = zeros(ndof,1);
    masterCol(masterList) = 1:nMaster;

    Ti = zeros(ndof,1); Tj = zeros(ndof,1); Tv = zeros(ndof,1); cnt = 0;
    for d = 1:ndof
        if isFixed(d)
            % up(d) already 0
        elseif isMaster(d)
            cnt=cnt+1; Ti(cnt)=d; Tj(cnt)=masterCol(d); Tv(cnt)=1;
        elseif isSlave(d)
            md = slaveMasterDof(d);
            if md ~= 0
                cnt=cnt+1; Ti(cnt)=d; Tj(cnt)=masterCol(md); Tv(cnt)=1;
            end
        end
    end
    T = sparse(Ti(1:cnt), Tj(1:cnt), Tv(1:cnt), ndof, nMaster);
end