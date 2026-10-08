function [V, dV, d2V] = legendre_basis(x, deg)
% Legendre polynomials P_0..P_deg (columns) at the points x in [-1,1] and their first and second derivatives
x = x(:);  V = zeros(numel(x), deg + 1);  dV = V;  d2V = V;
V(:, 1) = 1;
if deg >= 1, V(:, 2) = x;  dV(:, 2) = 1; end
for n = 2:deg, V(:, n+1) = ((2*n - 1)*x.*V(:, n) - (n - 1)*V(:, n-1))/n; end
for n = 1:deg-1
    dV(:, n+2)  = (2*n + 1)*V(:, n+1) + dV(:, n);
    d2V(:, n+2) = (2*n + 1)*dV(:, n+1) + d2V(:, n);
end
end
