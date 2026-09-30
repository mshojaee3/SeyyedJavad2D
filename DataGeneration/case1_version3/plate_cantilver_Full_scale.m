% =========================================================================
% 2D Plane Strain FEM: Periodic Inclusions Plate
% Domain: 1.0 m x 0.5 m, 6x3 Cells with circular inclusions
% Boundary: Left edge clamped (u=v=0), Right edge uniform shear traction
% Mesh: Conforming T6 quadratic elements (matches dataGeneration.m)
% =========================================================================

clear; close all; clc;

%% 1. Geometric & Material Parameters
Lx_tot = 1.0;          % Total length (m)
Ly_tot = 0.5;          % Total height (m)
nx = 12;                % Cells in X
ny = 6;                % Cells in Y

% Cell dimensions
Lx_c = Lx_tot / nx;    % Width of one cell: 1/6 m
Ly_c = Ly_tot / ny;    % Height of one cell: 1/6 m (square cells)
L_cell = min(Lx_c, Ly_c);

% Parameters specified
Rfrac    = 6/19;       % Inclusion radius / cell size
R        = Rfrac * L_cell;
meshFrac = 0.05;       % Element size / cell size
Hmax     = meshFrac * L_cell;

% Material properties (Plane Strain)
E_m = 70000.0; nu_m = 0.33;    % Matrix (e.g. MPa or N/mm^2)
E_i = 3500.0;  nu_i = 0.33;    % Inclusion
THICKNESS = 1.0;

% Downward shear traction on right edge (Force per unit length, e.g., N/m or N/mm)
Traction_Y = -1.0;    

fprintf('Building geometry and conforming mesh for %dx%d cells...\n', nx, ny);

%% 2. Geometry Construction with PDE Toolbox (Conforming Inclusions)
% Outer boundary rectangle:
% Code 3 (polygon), 4 vertices, x1..x4, y1..y4
gd = [3; 4; 0; Lx_tot; Lx_tot; 0; 0; 0; Ly_tot; Ly_tot];

names = {'R1'};
sf = 'R1';

nInclusions = nx * ny;
centers = zeros(nInclusions, 2);
cIdx = 0;

for j = 1:ny
    for i = 1:nx
        cIdx = cIdx + 1;
        cx = (i - 0.5) * Lx_c;
        cy = (j - 0.5) * Ly_c;
        centers(cIdx, :) = [cx, cy];
        
        % Circle: Code 1, center (cx, cy), radius R, padded to 10 rows
        gd = [gd, [1; cx; cy; R; 0; 0; 0; 0; 0; 0]];
        names{end+1} = sprintf('C%d', cIdx);
        sf = [sf '+' names{end}];
    end
end
sf = ['(' sf ')*R1'];

% CRUCIAL FIX: char(names)' converts cell array of strings into valid char matrix
[dl, ~] = decsg(gd, sf, char(names)');

% Build PDE model & quadratic triangular (T6) mesh
model = createpde(1);
geometryFromEdges(model, dl);
meshObj = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', 'quadratic');

nodes = meshObj.Nodes';          % (nNodes x 2)
elements = meshObj.Elements';    % (nElem x 6)
nNodes = size(nodes, 1);
nElem  = size(elements, 1);
nDOF   = 2 * nNodes;

fprintf('Mesh successfully generated: %d Nodes, %d Quadratic Elements (T6).\n', nNodes, nElem);


%% 3. Material Assignment (Matrix vs. Inclusion)
% Test element centroids to check if they lie inside any inclusion
elemCentroids = (nodes(elements(:,1), :) + nodes(elements(:,2), :) + nodes(elements(:,3), :)) / 3;
isIncElem = false(nElem, 1);
for k = 1:nInclusions
    distSq = (elemCentroids(:,1) - centers(k,1)).^2 + (elemCentroids(:,2) - centers(k,2)).^2;
    isIncElem = isIncElem | (distSq <= (R + 1e-7)^2);
end

% Plane strain constitutive matrices
Dm = planeStrainD(E_m, nu_m);
Di = planeStrainD(E_i, nu_i);

%% 4. Assembly of Global Stiffness Matrix (K)
% 3-point Gauss quadrature on unit triangle
r_gp = [2/3, 1/6, 1/6];
s_gp = [1/6, 2/3, 1/6];
w_gp = [1/6, 1/6, 1/6];

% Sparse triplet storage allocation (12x12 = 144 entries per T6 element)
iK = zeros(nElem * 144, 1);
jK = zeros(nElem * 144, 1);
sK = zeros(nElem * 144, 1);
entryCount = 0;

for e = 1:nElem
    eNodes = elements(e, :);
    xe = nodes(eNodes, 1);
    ye = nodes(eNodes, 2);
    
    % Pick constitutive law
    if isIncElem(e)
        D = Di;
    else
        D = Dm;
    end
    
    ke = zeros(12, 12);
    for g = 1:3
        r = r_gp(g);
        s = s_gp(g);
        t = 1 - r - s;
        
        % Quadratic T6 shape function derivatives w.r.t area coords (r, s)
        % N1=r(2r-1), N2=s(2s-1), N3=t(2t-1), N4=4rs, N5=4st, N6=4rt
        dNr = [4*r - 1;       0; -(4*t - 1);  4*s; -4*s;  4*(t - r)];
        dNs = [      0; 4*s - 1; -(4*t - 1);  4*r;  4*(t - s); -4*r];
        
        % Jacobian: [dx/dr, dy/dr; dx/ds, dy/ds]
        J = [dNr'; dNs'] * [xe, ye];
        detJ = det(J);
        invJ = inv(J);
        
        % Derivatives w.r.t physical coordinates (x, y)
        dNxy = invJ * [dNr'; dNs']; % (2 x 6)
        
        % Strain-displacement matrix B (3 x 12)
        B = zeros(3, 12);
        for a = 1:6
            B(1, 2*a-1) = dNxy(1, a);
            B(2, 2*a)   = dNxy(2, a);
            B(3, 2*a-1) = dNxy(2, a);
            B(3, 2*a)   = dNxy(1, a);
        end
        
        ke = ke + (B' * D * B) * detJ * w_gp(g) * THICKNESS;
    end
    
    % DOF indexing for current element
    eDOFs = reshape([2*eNodes-1; 2*eNodes], 1, 12);
    [Igrid, Jgrid] = meshgrid(eDOFs, eDOFs);
    idxRange = entryCount + (1:144);
    iK(idxRange) = Igrid(:);
    jK(idxRange) = Jgrid(:);
    sK(idxRange) = ke(:);
    entryCount = entryCount + 144;
end

K = sparse(iK, jK, sK, nDOF, nDOF);

%% 5. Boundary Conditions and Consistent Traction Loading
F = zeros(nDOF, 1);
TOL = 1e-6;

% 5.1 Left Edge (x = 0): Fully Clamped (u = 0, v = 0)
leftNodeIdx = find(abs(nodes(:, 1)) < TOL);
fixedDOFs = [2*leftNodeIdx - 1; 2*leftNodeIdx];
fixedDOFs = unique(fixedDOFs);

% 5.2 Right Edge (x = Lx_tot): Downward Traction loading
% Locate T6 edges that lie on the right boundary
rightEdges = [];
for e = 1:nElem
    en = elements(e, :);
    % T6 edges: 1-4-2, 2-5-3, 3-6-1
    t6_edges = [en(1), en(4), en(2);
                en(2), en(5), en(3);
                en(3), en(6), en(1)];
    for ed = 1:3
        edgeNodes = t6_edges(ed, :);
        if all(abs(nodes(edgeNodes, 1) - Lx_tot) < TOL)
            rightEdges = [rightEdges; edgeNodes]; %#ok<AGROW>
        end
    end
end

% Consistent 3-node Simpson line load: Int(N_i * Traction * ds)
% Shape functions for 1D quadratic edge: [-1, 1]
% Simpson weights for length L_e: [L_e/6, 4*L_e/6, L_e/6]
for ed = 1:size(rightEdges, 1)
    n1 = rightEdges(ed, 1);
    nm = rightEdges(ed, 2);
    n2 = rightEdges(ed, 3);
    
    Le = abs(nodes(n2, 2) - nodes(n1, 2));
    f_end = (Le / 6) * Traction_Y * THICKNESS;
    f_mid = (4 * Le / 6) * Traction_Y * THICKNESS;
    
    F(2*n1) = F(2*n1) + f_end;
    F(2*nm) = F(2*nm) + f_mid;
    F(2*n2) = F(2*n2) + f_end;
end

%% 6. Solution of Linear System
freeDOFs = setdiff(1:nDOF, fixedDOFs);
U = zeros(nDOF, 1);
U(freeDOFs) = K(freeDOFs, freeDOFs) \ F(freeDOFs);

u_x = U(1:2:end);
u_y = U(2:2:end);

%% 7. Midpoint Displacement of the Right Edge
rightNodes = find(abs(nodes(:, 1) - Lx_tot) < TOL);
[~, midSubIdx] = min(abs(nodes(rightNodes, 2) - (Ly_tot / 2)));
midNode = rightNodes(midSubIdx);

fprintf('\n===================================================\n');
fprintf('  RIGHT EDGE MIDPOINT RESULTS (x = %.4f, y = %.4f)\n', nodes(midNode,1), nodes(midNode,2));
fprintf('---------------------------------------------------\n');
fprintf('  Horizontal Displacement (u_x) = %+.6e m\n', u_x(midNode));
fprintf('  Vertical Displacement   (u_y) = %+.6e m\n', u_y(midNode));
fprintf('===================================================\n\n');

%% 8. Post-Processing: Mesh & Deformed Plot
scaleFactor = 100; % Visualization amplification
defNodes = nodes + scaleFactor * [u_x, u_y];

figure('Color', 'white', 'Position', [100, 100, 1100, 520]);

% Subplot 1: Microstructure Material Map
subplot(1, 2, 1); hold on; box on; axis equal;
title('Microstructure (Matrix & Circular Inclusions)', 'FontSize', 11);
% Draw matrix elements
patch('Faces', elements(~isIncElem, 1:3), 'Vertices', nodes, ...
      'FaceColor', [0.85 0.90 0.95], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.2);
% Draw inclusion elements
patch('Faces', elements(isIncElem, 1:3), 'Vertices', nodes, ...
      'FaceColor', [0.95 0.60 0.50], 'EdgeColor', [0.8 0.3 0.2], 'LineWidth', 0.2);
xlabel('x (m)'); ylabel('y (m)');
xlim([-0.05, Lx_tot + 0.05]); ylim([-0.05, Ly_tot + 0.05]);

% Subplot 2: Deformed Mesh with Vertical Displacement Contours
subplot(1, 2, 2); hold on; box on; axis equal;
title(sprintf('Deformed Configuration (Scale: %dx)', scaleFactor), 'FontSize', 11);

% Contour of u_y on deformed geometry
patch('Faces', elements(:, 1:3), 'Vertices', defNodes, ...
      'FaceVertexCData', u_y, 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(jet(256));
cb = colorbar;
ylabel(cb, 'Vertical Displacement v (m)', 'FontSize', 10);

% Mark the original and deformed midpoint
plot(nodes(midNode, 1), nodes(midNode, 2), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
plot(defNodes(midNode, 1), defNodes(midNode, 2), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

legend('', 'Original Midpoint', 'Deformed Midpoint', 'Location', 'southoutside');
xlabel('x (m)'); ylabel('y (m)');
xlim([-0.05, Lx_tot + 0.1]); ylim([-0.2, Ly_tot + 0.05]);

%% Helper Function
function D = planeStrainD(E, nu)
    c = E / ((1 + nu) * (1 - 2*nu));
    D = c * [ 1 - nu,     nu,            0;
                nu,   1 - nu,            0;
                 0,        0,   (1 - 2*nu)/2 ];
end
