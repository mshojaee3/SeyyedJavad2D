function M = micromorphic_solution(P, par)
% =========================================================================
%  MICROMORPHIC THEORY SOLUTION OF THE TEST PROBLEM (forward finite elements with the identified parameters)   M = micromorphic_solution(P, par)
% =========================================================================
%  Plate P.test, clamped at x = 0, traction t on x = Lx, for every N in NList (model and element: see micromorphic_fe.m):
%      gamma = eps_macro - psi,   W = 1/2 psi'C_psi(N) psi + 1/2 gamma'C_gamma(N) gamma + psi'C_couple(N) gamma + 1/2 l^2 G |grad psi|^2,   l = ellFactor*Ly/N
%      sigma = C_gamma gamma + C_couple' psi          virtual work:  K U = F ,   K = sum_p theta_p(N) Kp{p} + l^2 G Kg  (C_psi, C_gamma, C_couple: 9 parameters)
%  u = v = 0 on x = 0; psi free (natural condition l^2 G grad(psi) n = 0).  N outside the identified range = extrapolation.
%
%  OUTPUT M.N(iN): N, twoU = F'U, ux, uy (mean edge displacement), react = [Rx; Ry; M] of the clamp (K U - F), minEig of [C_psi C_c; C_c' C_gamma],
%                  Pbin, Sbin, Ebin: psi, sigma and eps_macro averaged on bins (nBinY x nBinX x 3 x nLoads);   M.xc, M.yc: bin centres, M.dB: bin size
%                  raw: the FE solution itself: U = nodal unknowns [u v psi11 psi22 psi12] of every node (nDof x nLoads), Gauss-point data Xg, Wg,
%                       Psi, Sig, Eps (single; nElem x 6 points [x 3 components x nLoads]);   M.mesh: nodes and elements (T6) of the micromorphic mesh
%                  (there is no smoothing step: the micromorphic fields are used as they are)
% =========================================================================
T = P.test;  Lx = T.Lx;  Ly = T.Ly;  Ns = T.NList;  nL = size(T.trc, 3);  
nby = P.fe.nBinY;  nbx = round(nby*Lx/Ly);  assert(abs(nbx*Ly/nby - Lx) < 1e-9, 'nBinY*Lx/Ly must be an integer.');
S = micromorphic_fe(P, 'test', P.fe.nyMesh);  nd = size(S.F, 1);  M.mesh = struct('nodes', S.nodes, 'elems', S.elems);
fprintf('  micromorphic: plate %.3g x %.3g, N = %s, mesh %d x %d\n', Lx, Ly, mat2str(Ns), S.nx, S.ny);
M.xc = ((1:nbx) - 0.5)*Ly/nby;  M.yc = ((1:nby) - 0.5)*Ly/nby;  M.dB = Ly/nby;  M.names = T.loads;
for iN = 1:numel(Ns)
    N = Ns(iN);  ell = par.ellFactor*Ly/N;
    if N < par.Nrange(1) || N > par.Nrange(2), fprintf('   note: N = %d is outside the identified range [%g, %g] (extrapolation)\n', N, par.Nrange); end
    Cg = stiffness_at(par.c_gamma, N);  Cc = stiffness_at(par.c_couple, N);  Cp = stiffness_at(par.c_psi, N);  M6 = [Cp, Cc; Cc.', Cg];  me = min(eig((M6 + M6.')/2));
    if me < 0, fprintf('   WARNING N = %d: energy matrix not positive definite (min eig %.3g): the solution may be unreliable\n', N, me); end
    th = [Cg(1,1), Cg(1,2), Cg(3,3), Cc(1,1), Cc(1,2), Cc(3,3), Cp(1,1), Cp(1,2), Cp(3,3)];
    K = ell^2*par.G*S.Kg;
    for p = 1:9, K = K + th(p)*S.Kp{p}; end
    U = zeros(nd, nL);  U(S.free, :) = K(S.free, S.free) \ S.F(S.free, :);
    Rf = K*U - S.F;  Rx = Rf(5*S.left - 4, :);  Ry = Rf(5*S.left - 3, :);
    [Xg, Wg, Pg, Sg, Eg] = S.gauss(U, Cg, Cc);
    raw = struct('U', U, 'Xg', single(Xg), 'Wg', single(Wg), 'Psi', single(Pg), 'Sig', single(Sg), 'Eps', single(Eg));
    r = struct('N', N, 'twoU', sum(S.F.*U, 1), 'ux', S.edgeX.'*U, 'uy', S.edgeY.'*U, 'minEig', me, ...
               'react', [sum(Rx, 1); sum(Ry, 1); -sum((S.nodes(S.left, 2) - Ly/2).*Rx, 1)], 'Pbin', zeros(nby, nbx, 3, nL), 'Sbin', zeros(nby, nbx, 3, nL), ...
               'Ebin', zeros(nby, nbx, 3, nL), 'raw', raw);
    for l = 1:nL
        r.Pbin(:, :, :, l) = bin_mean(Xg, Wg, Pg(:, :, :, l), Lx, nbx, nby);  r.Sbin(:, :, :, l) = bin_mean(Xg, Wg, Sg(:, :, :, l), Lx, nbx, nby);
        r.Ebin(:, :, :, l) = bin_mean(Xg, Wg, Eg(:, :, :, l), Lx, nbx, nby);
    end
    M.N(iN) = r;
    fprintf('   N = %d: l = %.4f, u_y(edge) = %s, 2U = %s\n', N, ell, mat2str(r.uy, 5), mat2str(r.twoU, 5));
end
end

function C = stiffness_at(c, N)
th = c(:, 1) + c(:, 2)/N + c(:, 3)/N^2;
C = [th(1) th(2) 0; th(2) th(1) 0; 0 0 th(3)];
end

function Z = bin_mean(Xg, Wg, vals, Lx, nbx, nby)
% area-weighted mean of Gauss-point values in the bins of a regular grid (square bins of size Lx/nbx)
d = Lx/nbx;  ix = min(max(floor(Xg(:,:,1)/d) + 1, 1), nbx);  iy = min(max(floor(Xg(:,:,2)/d) + 1, 1), nby);
idx = [iy(:), ix(:)];  W = accumarray(idx, Wg(:), [nby, nbx]);  Z = zeros(nby, nbx, 3);
for k = 1:3, v = vals(:, :, k);  Z(:, :, k) = accumarray(idx, Wg(:).*v(:), [nby, nbx])./max(W, eps); end
end
