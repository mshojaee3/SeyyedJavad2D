function Z = legendre_grid(C, nbx, nby)
% Smooth Legendre fields -> values at the centres of the nby x nbx bins of a regular grid over the plate (same convention as legendre_fit
% in fullscale_solution.m).   C: (degY+1) x (degX+1) x 3 x nLoads   ->   Z: nby x nbx x 3 x nLoads  (row 1 = bottom bin, y = 0)
s = ((1:nbx)' - 0.5)/nbx*2 - 1;  t = ((1:nby)' - 0.5)/nby*2 - 1;
Px = legendre_basis(s, size(C, 2) - 1);  Py = legendre_basis(t, size(C, 1) - 1);
sz = size(C);  nC = prod(sz(3:end));
C3 = reshape(C, sz(1), sz(2), nC);  Z3 = zeros(nby, nbx, nC);
for k = 1:nC, Z3(:, :, k) = Py*C3(:, :, k)*Px.'; end
Z = reshape(Z3, [nby, nbx, sz(3:end)]);
end
