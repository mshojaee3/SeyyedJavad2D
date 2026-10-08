function [nodes, elements] = mesh_structured(Nx, Ny, Lx, Ly)
% Nx x Ny squares, each split into two quadratic triangles (T6; node order: 3 corners, then mid-sides 12, 23, 31)
xs = linspace(0, Lx, 2*Nx + 1);  ys = linspace(0, Ly, 2*Ny + 1);
[X, Y] = meshgrid(xs, ys);
nodes = [reshape(X', [], 1), reshape(Y', [], 1)];
idx = @(i, j) j*(2*Nx + 1) + i + 1;
[a, b] = meshgrid(0:Nx-1, 0:Ny-1);
i = 2*a(:);  j = 2*b(:);
elements = [idx(i,j), idx(i+2,j),   idx(i+2,j+2), idx(i+1,j),   idx(i+2,j+1), idx(i+1,j+1); ...
            idx(i,j), idx(i+2,j+2), idx(i,j+2),   idx(i+1,j+1), idx(i+1,j+2), idx(i,j+1)];
end
