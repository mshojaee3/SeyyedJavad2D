function S = micromorphic_fe(P, which, ny)
% =========================================================================
%  MICROMORPHIC FINITE ELEMENT MODEL (assembled ONCE; shared by micromorphic_solution and the refinement in identification)
%  S = micromorphic_fe(P, 'test', ny)    plate P.test   |   S = micromorphic_fe(P, 'ident', ny)   plate P.ident,   ny = elements along the height
% =========================================================================
%  Quadratic triangles (T6), 5 dofs per node [u v psi11 psi22 psi12], clamped at x = 0, traction on x = Lx.
%      STANDARD MINDLIN FORM (first term built on the MACRO strain):
%      W = 1/2 eps'C_hom eps + 1/2 gamma'C_gamma gamma + eps'C_couple gamma + 1/2 l^2 G |grad psi|^2 ,   eps = sym grad u,   gamma = eps - psi
%  The energy is LINEAR in the parameters, so the stiffness matrix is  K = Khom + sum_p theta_p Kp{p} + l^2 G Kg  with  C = a B1 + b B2 + d B3  for each matrix:
%      Khom : fixed C_hom term (full C_hom matrix),   Kp{1..3} : a, b, d of C_gamma,   Kp{4..6} : a, b, d of C_couple,   Kg : gradient term (l = G = 1)
%  S.F (loads), S.free, S.left (clamped nodes), S.edgeX, S.edgeY (functionals: mean ux, uy of the edge x = Lx), S.nodes, S.elems,
%  S.gauss(U, Cg, Cc) -> [Xg, Wg, Pg, Sg, Eg]: Gauss-point positions, areas, psi, sigma = dW/deps = C_hom eps + C_gamma gamma + C_couple gamma + C_couple' eps,
%                        and the macro strain eps = sym grad u   (dimension 4 of Pg, Sg, Eg: load case)
%  Unknown ordering of U (nd x nLoads):  U(1:5:end) = u,  U(2:5:end) = v,  U(3:5:end) = psi11,  U(4:5:end) = psi22,  U(5:5:end) = psi12 (node by node)
% =========================================================================
T = P.(which);  Lx = T.Lx;  Ly = T.Ly;  trc = T.trc;  nL = size(trc, 3);  THK = P.material.thickness;
nx = round(ny*Lx/Ly);  assert(abs(nx*Ly/ny - Lx) < 1e-9, 'ny*Lx/Ly must be an integer (ny = %d).', ny);
[nodes, elems] = mesh_structured(nx, ny, Lx, Ly);  nd = 5*size(nodes, 1);

%% loads, clamp, edge functionals
TOL = 1e-6;  pairs = [1 4 2; 2 5 3; 3 6 1];  re = zeros(0, 3);
for ed = 1:3
    n1 = elems(:, pairs(ed,1));  nm = elems(:, pairs(ed,2));  n2 = elems(:, pairs(ed,3));
    m = abs(nodes(n1,1) - Lx) < TOL & abs(nodes(nm,1) - Lx) < TOL & abs(nodes(n2,1) - Lx) < TOL;
    re = [re; n1(m), nm(m), n2(m)]; %#ok<AGROW>
end
Le = abs(nodes(re(:,3), 2) - nodes(re(:,1), 2));  gp = [-sqrt(3/5), 0, sqrt(3/5)];  gw = [5/9, 8/9, 5/9];  F = zeros(nd, nL);
for q = 1:size(re, 1)
    nn = re(q, :).';  yy = nodes(nn, 2);
    for ig = 1:3
        xi = gp(ig);  Nq = [xi*(xi-1)/2, 1 - xi^2, xi*(xi+1)/2];  dN = [xi - 0.5, -2*xi, xi + 0.5];
        y = Nq*yy;  jac = abs(dN*yy);  s = 2*y/Ly - 1;  Pl = [1, s, (3*s^2 - 1)/2];
        for l = 1:nL
            F(5*nn - 4, l) = F(5*nn - 4, l) + Nq.'*(Pl*trc(1, :, l).')*jac*gw(ig)*THK;
            F(5*nn - 3, l) = F(5*nn - 3, l) + Nq.'*(Pl*trc(2, :, l).')*jac*gw(ig)*THK;
        end
    end
end
left = find(abs(nodes(:, 1)) < TOL);
S.edgeX = zeros(nd, 1);  S.edgeY = zeros(nd, 1);                          % mean over the edge x = Lx (Simpson on the element sides)
for q = 1:size(re, 1)
    w = Le(q)/6/sum(Le)*[1 4 1];
    S.edgeX(5*re(q, :) - 4) = S.edgeX(5*re(q, :) - 4) + w.';  S.edgeY(5*re(q, :) - 3) = S.edgeY(5*re(q, :) - 3) + w.';
end
S.F = F;  S.free = setdiff((1:nd).', [5*left - 4; 5*left - 3]);  S.left = left;  S.nodes = nodes;  S.elems = elems;  S.nx = nx;  S.ny = ny;

%% stiffness matrices  K_q = int L' Q_q L dOmega,  L = [BU; BU - NP; BG]  maps the element dofs to [eps; gamma; grad psi]
%  (eps = macro strain, gamma = eps - psi).  Matrices: 1..3 C_gamma (a,b,d), 4..6 C_couple (a,b,d), 7 gradient term, 8 fixed C_hom term
Chom = P.material.Chom;  Chom = (Chom + Chom.')/2;
B = {diag([1 1 0]), [0 1 0; 1 0 0; 0 0 0], diag([0 0 1])};  nQ = 8;  Qs = cell(1, nQ);
for q = 1:nQ, Qs{q} = zeros(12); end
for p = 1:3
    Qs{p}(4:6, 4:6) = B{p};                                                % C_gamma : 1/2 gamma' B gamma
    Qs{3+p}(1:3, 4:6) = B{p};  Qs{3+p}(4:6, 1:3) = B{p}.';                 % C_couple: eps' B gamma
end
Qs{7}(7:12, 7:12) = eye(6);                                                % gradient term (l = G = 1)
Qs{8}(1:3, 1:3) = Chom;                                                    % 1/2 eps' C_hom eps (fixed)
[gr, gs, gwt] = gauss6();  nE = size(elems, 1);  K = cell(1, nQ);
for q = 1:nQ, K{q} = sparse(nd, nd); end
for c0 = 1:400:nE                                                          % chunks of elements (memory)
    ee = c0:min(c0 + 399, nE);  nC = numel(ee);  iK = zeros(900, nC);  jK = iK;  vK = cell(1, nQ);
    for q = 1:nQ, vK{q} = zeros(900, nC); end
    for k = 1:nC
        en = elems(ee(k), :);  Ke = zeros(30, 30, nQ);
        for g = 1:6
            [BU, NP, BG, detJ] = gauss_matrices(gr(g), gs(g), nodes(en, 1), nodes(en, 2));
            L = [BU; BU - NP; BG];
            for q = 1:nQ, Ke(:, :, q) = Ke(:, :, q) + (L.'*Qs{q}*L)*(detJ*gwt(g)); end
        end
        dofs = reshape((5*(en(:) - 1) + (1:5)).', 1, []);  [jj, ii] = meshgrid(dofs, dofs);
        iK(:, k) = ii(:);  jK(:, k) = jj(:);
        for q = 1:nQ, v = Ke(:, :, q);  vK{q}(:, k) = v(:); end
    end
    for q = 1:nQ, K{q} = K{q} + sparse(iK(:), jK(:), vK{q}(:), nd, nd); end
end
S.Kp = K(1:6);  S.Kg = K{7};  S.Khom = K{8};
S.gauss = @(U, Cg, Cc) gauss_fields(nodes, elems, U, Chom, Cg, Cc);
end

function [BU, NP, BG, detJ, N] = gauss_matrices(r, s, xe, ye)
% matrices at one Gauss point:  BU*ue = eps_macro,  NP*ue = psi,  BG*ue = [psi,x; psi,y]
t = 1 - r - s;
N   = [r*(2*r-1); s*(2*s-1); t*(2*t-1); 4*r*s; 4*s*t; 4*t*r];
dNr = [4*r - 1; 0; -(4*t - 1); 4*s; -4*s; 4*(t - r)];
dNs = [0; 4*s - 1; -(4*t - 1); 4*r; 4*(t - s); -4*r];
J = [dNr.'*xe, dNr.'*ye; dNs.'*xe, dNs.'*ye];  detJ = det(J);  dN = J \ [dNr.'; dNs.'];
BU = zeros(3, 30);  NP = zeros(3, 30);  BG = zeros(6, 30);
for a = 1:6
    c0 = 5*(a - 1);
    BU(:, c0 + (1:2)) = [dN(1,a) 0; 0 dN(2,a); dN(2,a) dN(1,a)];
    NP(:, c0 + (3:5)) = N(a)*eye(3);
    BG(1:3, c0 + (3:5)) = dN(1,a)*eye(3);  BG(4:6, c0 + (3:5)) = dN(2,a)*eye(3);
end
end

function [r, s, w] = gauss6()
% 6-point rule of degree 4 on the reference triangle
a = 0.445948490915965;  b = 0.091576213509771;  wa = 0.223381589678011;  wb = 0.109951743655322;
r = [a, a, 1 - 2*a, b, b, 1 - 2*b];  s = [a, 1 - 2*a, a, b, 1 - 2*b, b];  w = 0.5*[wa wa wa wb wb wb];
end

function [Xg, Wg, Pg, Sg, Eg] = gauss_fields(nodes, elems, U, Chom, Cg, Cc)
% psi, macro strain eps and sigma = dW/deps = C_hom eps + C_gamma gamma + C_couple gamma + C_couple' eps at the 6 Gauss points of every element (4th dimension: load case)
nE = size(elems, 1);  nL = size(U, 2);  [gr, gs, gw] = gauss6();
Xg = zeros(nE, 6, 2);  Wg = zeros(nE, 6);  Pg = zeros(nE, 6, 3, nL);  Sg = Pg;  Eg = Pg;
for e = 1:nE
    en = elems(e, :);  xe = nodes(en, 1);  ye = nodes(en, 2);  ue = U(reshape((5*(en(:) - 1) + (1:5)).', 1, []), :);
    for g = 1:6
        [BU, NP, ~, detJ, N] = gauss_matrices(gr(g), gs(g), xe, ye);
        Xg(e, g, 1) = N.'*xe;  Xg(e, g, 2) = N.'*ye;  Wg(e, g) = detJ*gw(g);
        epsm = BU*ue;  psi = NP*ue;  gam = epsm - psi;  sig = Chom*epsm + Cg*gam + Cc*gam + Cc.'*epsm;
        for l = 1:nL, Pg(e, g, :, l) = psi(:, l);  Sg(e, g, :, l) = sig(:, l);  Eg(e, g, :, l) = epsm(:, l); end
    end
end
end
