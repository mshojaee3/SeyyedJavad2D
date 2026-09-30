%% fast_micromorphic_5dof_cantilever.m
% =========================================================================
% Direct 5-DOF Micromorphic Finite Element Solver for Cantilever Plate
% DOFs   : [u, v, psi11, psi22, psi12] per node
% Plate  : 1.0 m x 0.5 m, Clamped at Left (x=0), Downward traction at Right (x=1)
% =========================================================================
clear; clc; close all;

tic;
fprintf('=== [1] Setting Pre-computed Inputs ===\n');

%% 1. DIRECT INPUTS
% -------------------------------------------------------------------------
Chom =  [4.7854, 1.7582, 0;
        1.7582, 4.7875, 0;
        0,       0,     1.0305] * 1.0e4;

alpha_m = 1.0000;
alpha_h =  1.26929 ;
alpha_c =  1.06446 ;
alpha_l = 0.0959719;

ell  = 1.0 / 2.0;
ell2 = ell^2;

fprintf('  Chom provided (3x3).\n');
fprintf('  Alphas: alpha_m=%.4f, alpha_h=%.4f, alpha_c=%.4f, alpha_l=%.4f\n', ...
    alpha_m, alpha_h, alpha_c, alpha_l);
fprintf('  Characteristic length ell = %.4f m\n', ell);

%% 2. PLATE GEOMETRY, MESH AND DATA EXTRACTION
% -------------------------------------------------------------------------
fprintf('=== [2] Creating Geometry and Mesh ===\n');

Lx = 1.0;          % Length (m)
Ly = 0.5;          % Height (m)
traction_y = -1.0; % Downward traction (N/m)

model = createpde();

% Decomposed geometry matrix for rectangle [0, Lx] x [0, Ly]
% Row 1: 3 (decomposed geometry format)
% Row 2: 4 (number of edges)
% Rows 3-6: x-coordinates of the four corners
% Rows 7-10: y-coordinates of the four corners
rect = [3; 4; 0; Lx; Lx; 0; 0; 0; Ly; Ly];
gd = decsg(rect);                 % decompose into edge segments
geometryFromEdges(model, gd);     % build geometry from edge matrix

% Generate quadratic triangular mesh
msh = generateMesh(model, 'Hmax', 0.04, 'GeometricOrder', 'quadratic');

% ---- Extract mesh data ----
nodes    = msh.Nodes';       % nNode x 2
elements = msh.Elements';    % nElem x 6  (T6 quadratic triangles)
nNode    = size(nodes, 1);
nElem    = size(elements, 1);

fprintf('  Mesh: %d nodes, %d quadratic triangular elements\n', nNode, nElem);

%% 3. ELEMENT STIFFNESS MATRIX & SYSTEM ASSEMBLY
% -------------------------------------------------------------------------
fprintf('=== [3] Assembling 5-DOF Coupled Global Stiffness Matrix ===\n');

nDOF_total = 5 * nNode;
K = sparse(nDOF_total, nDOF_total);
F = zeros(nDOF_total, 1);

% 6-point Gaussian quadrature for quadratic triangle
gp_r = [1/6, 2/3, 1/6, 1/2, 0.0, 1/2];
gp_s = [1/6, 1/6, 2/3, 1/2, 1/2, 0.0];
gp_w = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6] * 0.5;

% Metric for curvature energy
A6 = blkdiag(Chom, Chom);

% Constitutive coupling matrices
C_uu   = alpha_h * Chom;
C_up   = (alpha_c - alpha_h) * Chom;
C_pp   = (alpha_m + alpha_h - 2 * alpha_c) * Chom;
C_grad = alpha_l * ell2 * A6;

for e = 1:nElem
    elemNodes = elements(e, :);
    xe = nodes(elemNodes, 1);
    ye = nodes(elemNodes, 2);

    Ke = zeros(30, 30); % 6 nodes * 5 DOFs

    for gp = 1:6
        r = gp_r(gp);
        s = gp_s(gp);
        w = gp_w(gp);
        t = 1 - r - s;

        % T6 quadratic shape functions
        N = [t*(2*t - 1), r*(2*r - 1), s*(2*s - 1), 4*r*t, 4*r*s, 4*s*t];

        % Derivatives wrt reference coordinates
        dN_dr = [1 - 4*t, 4*r - 1, 0,       4*(t - r), 4*s,       -4*s];
        dN_ds = [1 - 4*t, 0,       4*s - 1, -4*r,      4*r,       4*(t - s)];

        J = [dN_dr * xe, dN_dr * ye;
             dN_ds * xe, dN_ds * ye];
        detJ = det(J);
        invJ = inv(J);

        % Cartesian derivatives [2 x 6]
        dN_dxy = invJ * [dN_dr; dN_ds];

        % 1. Macro strain-displacement Bu [3 x 12]
        Bu = zeros(3, 12);
        for a = 1:6
            Bu(:, (2*a-1):(2*a)) = [dN_dxy(1,a), 0;
                                    0,           dN_dxy(2,a);
                                    dN_dxy(2,a), dN_dxy(1,a)];
        end

        % 2. Micro-strain interpolation Npsi [3 x 18]
        Npsi = zeros(3, 18);
        for a = 1:6
            Npsi(:, (3*a-2):(3*a)) = diag([1, 1, 2]) * N(a);
        end

        % 3. Micro-strain gradient Bgrad [6 x 18]
        Bgrad = zeros(6, 18);
        for a = 1:6
            dNx = dN_dxy(1, a);
            dNy = dN_dxy(2, a);
            Bgrad(1:3, (3*a-2):(3*a)) = diag([dNx, dNx, 2*dNx]);
            Bgrad(4:6, (3*a-2):(3*a)) = diag([dNy, dNy, 2*dNy]);
        end

        % Coupled stiffness sub-matrices
        Ke_uu = Bu'   * C_uu   * Bu;
        Ke_up = Bu'   * C_up   * Npsi;
        Ke_pp = Npsi' * C_pp   * Npsi + Bgrad' * C_grad * Bgrad;

        Ke_local = [Ke_uu,    Ke_up;
                    Ke_up',   Ke_pp];

        Ke = Ke + Ke_local * (detJ * w);
    end

    % Global DOF mapping for element: [u, v] then [psi11, psi22, psi12]
    dof_u   = reshape([2*elemNodes-1; 2*elemNodes], 1, []);
    dof_psi = reshape([2*nNode + 3*elemNodes-2;
                       2*nNode + 3*elemNodes-1;
                       2*nNode + 3*elemNodes], 1, []);
    elem_dofs = [dof_u, dof_psi];

    K(elem_dofs, elem_dofs) = K(elem_dofs, elem_dofs) + Ke;
end

%% 4. BOUNDARY CONDITIONS & LOADS
% -------------------------------------------------------------------------
fprintf('=== [4] Applying Boundary Conditions and Loads ===\n');

tol = 1e-5;

% A. Transverse downward traction on right edge (x = Lx)
% Edge-node connectivity in T6 triangle: [1 4 2], [2 5 3], [3 6 1]
for e = 1:nElem
    en = elements(e, :);
    edges = [en([1 4 2]); en([2 5 3]); en([3 6 1])];
    for ed = 1:3
        edge_nodes = edges(ed, :);
        if all(abs(nodes(edge_nodes, 1) - Lx) < tol)
            L_edge = abs(nodes(edge_nodes(3), 2) - nodes(edge_nodes(1), 2));
            % Simpson weights on 3-node quadratic edge: (L/6) * [1; 4; 1]
            f_local = (traction_y * L_edge / 6) * [1; 4; 1];
            F(2 * edge_nodes) = F(2 * edge_nodes) + f_local;
        end
    end
end

% B. Fully clamped BCs on left edge (x = 0)
left_nodes = find(abs(nodes(:, 1)) < tol);
fixed_dofs = [
    2 * left_nodes - 1;          % u = 0
    2 * left_nodes;              % v = 0
    2*nNode + 3*left_nodes - 2;  % psi11 = 0
    2*nNode + 3*left_nodes - 1;  % psi22 = 0
    2*nNode + 3*left_nodes       % psi12 = 0
];
fixed_dofs = unique(fixed_dofs);

all_dofs  = 1:nDOF_total;
free_dofs = setdiff(all_dofs, fixed_dofs);

%% 5. SOLUTION & RESULTS EXTRACTION
% -------------------------------------------------------------------------
fprintf('=== [5] Solving Sparse Linear System ===\n');

U = zeros(nDOF_total, 1);
U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);

% Extract fields
ux = U(1:2:2*nNode);
uy = U(2:2:2*nNode);

% Locate midpoint of right edge: (x = Lx, y = Ly/2)
mid_pt = [Lx, Ly / 2];
dist = sum((nodes - mid_pt).^2, 2);
[~, mid_idx] = min(dist);

solve_time = toc;
fprintf('\n======================================================\n');
fprintf('                 FINAL RESULTS\n');
fprintf('======================================================\n');
fprintf('Solution completed in: %.2f seconds\n', solve_time);
fprintf('Identified Alphas:\n');
fprintf('  alpha_m = %.6f\n', alpha_m);
fprintf('  alpha_h = %.6f\n', alpha_h);
fprintf('  alpha_c = %.6f\n', alpha_c);
fprintf('  alpha_l = %.6f\n', alpha_l);
fprintf('------------------------------------------------------\n');
fprintf('Right edge midpoint node: #%d at (x = %.4f m, y = %.4f m)\n', ...
    mid_idx, nodes(mid_idx, 1), nodes(mid_idx, 2));
fprintf('  Horizontal Displacement (u_x) : %+.8e m\n', ux(mid_idx));
fprintf('  Vertical Deflection     (u_y) : %+.8e m\n', uy(mid_idx));
fprintf('======================================================\n');

%% 6. DEFLECTION CONTOUR PLOT
figure('Name', 'Micromorphic Homogenized Cantilever Deflection', 'Color', 'w');
pdeplot(model, 'XYData', uy, 'Mesh', 'off', 'ColorMap', 'jet');
title(sprintf('Vertical Deflection v(x,y) [\\alpha_m=%.2f, \\alpha_h=%.2f, \\alpha_c=%.2f, \\alpha_l=%.2f]', ...
    alpha_m, alpha_h, alpha_c, alpha_l));
xlabel('X (m)'); ylabel('Y (m)'); axis equal; colorbar;

% --- Plot of deformed shape (optional)
figure('Name', 'Deformed Shape (scaled)', 'Color', 'w');
scale = 1e2;   % exaggerate deflection for visualization
pdeplot(model, 'XYData', uy, 'Mesh', 'on', 'ColorMap', 'parula');
hold on;
% Overlay scaled deformed nodes
Xdef = nodes(:,1) + scale*ux;
Ydef = nodes(:,2) + scale*uy;
plot(Xdef, Ydef, 'k.', 'MarkerSize', 4);
title(sprintf('Deformed shape (x%.0f), v_{max}=%.3e m', scale, min(uy)));
axis equal; colorbar;