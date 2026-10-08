function par = identification(P, D)
% =========================================================================
%  IDENTIFICATION OF THE MICROMORPHIC PARAMETERS   par = identification(P, D)
%  D = fullscale_solution(P, 'ident')  (eps_macro of the homogeneous plate, psi and U_het of the heterogeneous plates)
% =========================================================================
%  MODEL (plane strain, Voigt vectors [11; 22; 12]):
%      gamma = eps_macro - psi
%      W = 1/2 psi'C_psi psi + 1/2 gamma'C_gamma gamma + psi'C_couple gamma + 1/2 l^2 G |grad psi|^2,   l = ellFactor*Ly/N
%      sigma = C_gamma gamma + C_couple' psi
%      C_psi(N), C_gamma(N), C_couple(N) = [a b 0; b a 0; 0 0 d]  with the SAME parameterization  theta(N) = c0 + c1/N + c2/N^2  for every entry;
%      each matrix is 'free' (identified) or 'hom' (= C_hom), see parameters.m.   Unknowns: 9 coefficients per free matrix + G.
%
%  STEP A - RESIDUALS OF THE THEORY (linear in the unknowns, weak form so that no derivative of a fitted field is needed), stacked for all N and loads:
%      macro  : div(sigma) = 0                         (test functions: bubbles in x and y)
%      micro  : dW/dpsi = C_psi psi - (C_gamma - C_couple) gamma - C_couple' psi - l^2 G Lap(psi) = 0
%      trac   : sigma n = t on x = Lx, sigma n = 0 on y = 0, Ly         (test functions: Legendre polynomials along the edge)
%      mu     : l^2 G (grad psi) n = 0 on the edges
%      energy : int W dOmega = U_het = 1/2 F'U of the centred heterogeneous plate
%      force  : int sigma_11 dy = Fx,  int sigma_12 dy = Fy      on the sections x = s Lx
%      moment : int (y - Ly/2) sigma_11 dy = -M_ext             on the sections x = s Lx
%  CONSTRAINTS (linear): 0 <= theta(N) <= boundFactor*C_hom,ij;  theta monotone in N (increasing/decreasing patterns, best kept);  theta(N -> inf) = C_hom,ij;
%      G >= 0;  mode by mode (sum, diff, shear) with A = C_psi, B = C_couple, C = C_gamma:  0 <= B <= A - margin, B <= C - margin  (positive definite energy)
%      and, for N >= Nhom,  B >= (1-homTol) C_hom,  A <= (1+homTol) C_hom  (long-wave stiffness B <= C_eff <= min(A,C) close to C_hom).
%  STEP B - REFINEMENT WITH THE FORWARD SOLUTION (model.refine): Levenberg-Marquardt on  J = ||A c - b||^2/||b||^2 + w*mean(r_u^2 + r_E^2),  r_u, r_E = relative
%      errors of the dominant edge displacement and of the energy 2U = F'U of the FORWARD micromorphic solution of the identification plate w.r.t. the
%      full-scale plate (all N and loads), same constraints, exact Jacobian (K is linear in the parameters, adjoint for the displacement).
%
%  OUTPUT par.c_psi, par.c_gamma, par.c_couple (3x3: rows a, b, d; columns c0, c1, c2; stiffness units), par.G, par.ellFactor, par.Nrange, par.minEig, par.relRes
% =========================================================================
fprintf('  identification: residual rows ...\n');
Id = P.ident;  mo = P.model;  rs = P.residual;  mat = P.material;  Chom = mat.Chom;  th0 = Chom(1, 1);
Ns = Id.NList;  nN = numel(Ns);  nL = size(Id.trc, 3);
hom3 = [(Chom(1,1) + Chom(2,2))/2; Chom(1,2); Chom(3,3)]/th0;        % a, b, d of C_hom in units of C_hom,11
ref9 = [hom3; hom3; hom3];  nS = 9;  nc = 3*nS + 1;                   % parameters q = 1..3 C_gamma, 4..6 C_couple, 7..9 C_psi (a, b, d); unknowns c = [c_q0 c_q1 c_q2 (q = 1..9), G]
xref = [ref9; mo.G0/th0];
fr = [strcmpi(mo.Cgamma, 'free'), strcmpi(mo.Ccouple, 'free'), strcmpi(mo.Cpsi, 'free')];       % which matrices are identified
freeQ = repelem(fr, 3);
Q = quad_setup(Id, rs);

%% A. residual rows of every N and load, in the joint unknowns
AJ = zeros(0, nc);  bJ = zeros(0, 1);  blocks = cell(nN, nL);
for iN = 1:nN
    N = Ns(iN);  phi = [1, 1/N, 1/N^2];
    for l = 1:nL
        blk = residual_rows(Id, rs, mo, Chom, D.hom.epsCoef(:,:,:,l), D.het(iN).psiCoef(:,:,:,l), Id.trc(:,:,l), D.het(iN).twoU(l)/2/mat.thickness, N, Q, xref);
        blocks{iN, l} = blk;
        for g = 1:numel(blk)
            A = blk(g).A;  Aj = zeros(size(A, 1), nc);
            for q = 1:nS, Aj(:, 3*(q-1) + (1:3)) = A(:, q)*phi; end
            Aj(:, nc) = A(:, nS + 1);
            AJ = [AJ; Aj];  bJ = [bJ; blk(g).b]; %#ok<AGROW>
        end
    end
end

%% constraints
Ngrid = unique([Ns(1):0.25:Ns(end), Ns]);  Ain = zeros(0, nc);  bin = zeros(0, 1);
for N = Ngrid
    for q = 1:nS
        r = zeros(1, nc);  r(3*(q-1) + (1:3)) = [1, 1/N, 1/N^2];
        Ain = [Ain; r; -r];  bin = [bin; mo.boundFactor*ref9(q); 0]; %#ok<AGROW>     % theta <= bound,  -theta <= 0
    end
end
if mo.stability                                                                       % stability + closeness to C_hom, mode by mode (all linear)
    Tm = [1 1 0; 1 -1 0; 0 0 1];  Am = Tm*hom3;                                       % modes sum = a+b, diff = a-b, shear = d and their C_hom values
    Nst = unique([Ngrid, Ns(end)*[1.5 2 3 4 6 10 20 50 100]]);
    for N = Nst
        phi = [1, 1/N, 1/N^2];
        for k = 1:3
            rowM = zeros(3, nc);                                                      % mode k of C_gamma (1), C_couple (2), C_psi (3) at N
            for m = 1:3, for q = 1:3, rowM(m, 3*((m-1)*3 + q - 1) + (1:3)) = Tm(k, q)*phi; end, end
            rowC = rowM(1, :);  rowB = rowM(2, :);  rowA = rowM(3, :);                % C = C_gamma, B = C_couple, A = C_psi
            mg = mo.margin*(N <= Ns(end));  mgBA = mg*(fr(3) || fr(2));  mgBC = mg*(fr(1) || fr(2));
            Ain = [Ain; -rowB; rowB - rowA; rowB - rowC]; %#ok<AGROW>                % B >= 0,  B <= A - margin,  B <= C - margin
            bin = [bin; 0; -mgBA*Am(k); -mgBC*Am(k)]; %#ok<AGROW>
            if N >= mo.Nhom
                Ain = [Ain; -rowB; rowA]; %#ok<AGROW>                                 % B >= (1-homTol) C_hom,  A <= (1+homTol) C_hom
                bin = [bin; -(1 - mo.homTol)*Am(k); (1 + mo.homTol)*Am(k)]; %#ok<AGROW>
            end
        end
    end
end
Aeq = zeros(0, nc);  beq = zeros(0, 1);                                               % limit:  c0 = C_hom  for every free matrix
if strcmpi(mo.limit, 'all')
    for q = find(freeQ), r = zeros(1, nc);  r(3*(q-1) + 1) = 1;  Aeq = [Aeq; r];  beq = [beq; ref9(q)]; end %#ok<AGROW>
end
lb = -inf(nc, 1);  lb(nc) = 0;  ub = inf(nc, 1);  ub(nc) = 10^mo.Grange(2);           % G >= 0 (and G/C_hom,11 <= 10^Grange(2) in the refinement)

% matrices fixed to C_hom ('hom') are not unknowns: c0 = C_hom, c1 = c2 = 0 -> eliminate them (known part moved to the right-hand sides)
colFree = true(nc, 1);  cfix = zeros(nc, 1);
for q = find(~freeQ), colFree(3*(q-1) + (1:3)) = false;  cfix(3*(q-1) + 1) = ref9(q); end
var = find(colFree);  fix = find(~colFree);
bJr = bJ - AJ(:, fix)*cfix(fix);  AJr = AJ(:, var);
binr = bin - Ain(:, fix)*cfix(fix);  Ainr = Ain(:, var);  keep = any(abs(Ainr) > 0, 2);  Ainr = Ainr(keep, :);  binr = binr(keep);
Aeqr = Aeq(:, var);  beqr = beq;  if isempty(Aeq), Aeqr = [];  beqr = []; end
lbr = lb(var);  ubr = ub(var);
expand = @(cv) accumulate_c(cv, var, cfix, nc);

%% B. constrained least squares.  Monotone: dtheta/dN = -(c1 N + 2 c2)/N^3 -> for all N >= Nmin:  s=+1 (decreasing): c1 >= 0, c1 Nmin + 2 c2 >= 0
fprintf('  identification: constrained least squares ...\n');
[Qm, Rm] = qr(AJr, 0);  dm = Qm.'*bJr;                                % ||AJ c - bJ||^2 = ||Rm c - dm||^2 + const  (reduced unknowns)
[cv, AmBest] = best_pattern(Rm, dm, Ainr, binr, Aeqr, beqr, lbr, Ns(1), nS, nc, find(freeQ), var, mo.monotone);

%% C. refinement with the forward displacement and energy
if mo.refine
    fprintf('  identification: refinement with the forward displacement and energy ...\n');
    cref = zeros(nc, 1);  cref(1:3:3*nS) = ref9;  cref(nc) = mo.G0/th0;   % homogeneous reference: theta = C_hom, G = G0
    scale = norm(abs(AJ)*cref);                                           % size of the terms of the residual equations (normalisation of the residual part of J)
    [cv, ~, rbest] = refine(P, D, AJr, bJr, cv, [Ainr; AmBest], [binr; zeros(size(AmBest, 1), 1)], Aeqr, beqr, lbr, ubr, var, cfix, nc, scale);
    fprintf('   forward errors (%%) of the identification plate, displacement | energy  (rows N, columns loads):\n');
    for iN = 1:nN
        idx = (iN-1)*nL + (1:nL);  fprintf('   N = %d:  %s  |  %s\n', Ns(iN), sprintf('%+7.2f', 100*rbest(idx)), sprintf('%+7.2f', 100*rbest(nN*nL + idx)));
    end
end
c = expand(cv);

%% D. parameters and report
cred = reshape(c(1:3*nS), 3, nS).';                                    % nS x 3, units of C_hom,11
par.c_gamma = cred(1:3, :)*th0;  par.c_couple = cred(4:6, :)*th0;  par.c_psi = cred(7:9, :)*th0;  par.G = c(nc)*th0;
par.ellFactor = mo.ellFactor;  par.Nrange = [Ns(1), Ns(end)];  par.relMisfit = norm(AJ*c - bJ)/norm(bJ);
fprintf('  relative misfit of all residual rows: %.4g\n', par.relMisfit);
nm = {'C_gamma a (C11=C22)', 'C_gamma b (C12)', 'C_gamma d (C33)', 'C_couple a (C11=C22)', 'C_couple b (C12)', 'C_couple d (C33)', 'C_psi a (C11=C22)', 'C_psi b (C12)', 'C_psi d (C33)'};
fprintf('  identified  theta(N) = c0 + c1/N + c2/N^2:\n');
for q = 1:nS, fprintf('   %-22s c0 = %10.1f  c1 = %11.1f  c2 = %11.1f   (C_hom = %.1f)%s\n', nm{q}, cred(q, :)*th0, ref9(q)*th0, ternary(freeQ(q), '', '  [fixed = C_hom]')); end
fprintf('   G = %.5g\n', par.G);
Nshow = unique([Ns, 10, 20, 50]);  minVal = inf;  Tm = [1 1 0; 1 -1 0; 0 0 1];  Am = Tm*hom3*th0;
fprintf('   N     C_psi:  a      b      d  |  C_gamma:  a      b      d  |  C_couple:  a      b      d  | min eig | long-wave C_eff/C_hom (sum, diff, shear)\n');
for N = Nshow
    Cg = stiffness_at(par.c_gamma, N);  Cc = stiffness_at(par.c_couple, N);  Cp = stiffness_at(par.c_psi, N);  M6 = [Cp, Cc; Cc.', Cg];
    v = [Cp(1,1) Cp(1,2) Cp(3,3) Cg(1,1) Cg(1,2) Cg(3,3) Cc(1,1) Cc(1,2) Cc(3,3)];  minVal = min([minVal, v]);
    Cm = Tm*[Cg(1,1); Cg(1,2); Cg(3,3)];  Bm = Tm*[Cc(1,1); Cc(1,2); Cc(3,3)];  Apm = Tm*[Cp(1,1); Cp(1,2); Cp(3,3)];  Ce = Cm - (Cm - Bm).^2./(Apm + Cm - 2*Bm);
    fprintf('   %-4d %8.0f %6.0f %6.0f  %9.0f %6.0f %6.0f  %9.0f %6.0f %6.0f  %8.3g   %6.3f %6.3f %6.3f\n', N, v, min(eig((M6 + M6.')/2)), Ce./Am);
end
fprintf('   all C_ij(N) >= 0: %s   (smallest value %.1f)\n', ternary(minVal >= -1e-6*th0, 'yes', 'NO'), minVal);
par.minEig = zeros(1, nN);
for iN = 1:nN
    Cg = stiffness_at(par.c_gamma, Ns(iN));  Cc = stiffness_at(par.c_couple, Ns(iN));  Cp = stiffness_at(par.c_psi, Ns(iN));  M6 = [Cp, Cc; Cc.', Cg];  par.minEig(iN) = min(eig((M6 + M6.')/2));
end
if any(par.minEig < 0), fprintf('   NOTE: the energy matrix is not positive definite for some N: the forward solution may be unreliable.\n'); end
groups = fieldnames(rs.weights).';  par.relRes = zeros(nN, numel(groups));  par.groups = groups;
for iN = 1:nN
    Cg = stiffness_at(par.c_gamma, Ns(iN));  Cc = stiffness_at(par.c_couple, Ns(iN));  Cp = stiffness_at(par.c_psi, Ns(iN));
    x = [Cg(1,1); Cg(1,2); Cg(3,3); Cc(1,1); Cc(1,2); Cc(3,3); Cp(1,1); Cp(1,2); Cp(3,3); par.G]/th0;
    for l = 1:nL
        for g = 1:numel(blocks{iN, l})
            bl = blocks{iN, l}(g);  k = find(strcmp(groups, bl.group));
            par.relRes(iN, k) = par.relRes(iN, k) + norm(bl.A*x - bl.b)/(norm(abs(bl.A)*abs(x)) + norm(bl.b))/nL;
        end
    end
end
fprintf('   relative residual per group and N (mean over loads):\n   %3s', 'N');  fprintf(' %9s', groups{:});  fprintf('\n');
for iN = 1:nN, fprintf('   %3d', Ns(iN));  fprintf(' %9.4f', par.relRes(iN, :));  fprintf('\n'); end
end

function [c, AmBest] = best_pattern(Rm, dm, Ain, bin, Aeq, beq, lb, Nmin, nS, nc, fq, var, monotone)
% constrained least squares  min ||Rm c - dm||^2  with the monotone direction s_q (+1 decreasing, -1 increasing) of every FREE parameter chosen by
% search: all patterns if there are at most 6 free parameters, otherwise a local search (flip one direction at a time, from three starting patterns)
opts = optimoptions('lsqlin', 'Display', 'off');  nF = numel(fq);
if ~monotone
    c = lsqlin(Rm, dm, Ain, bin, Aeq, beq, lb, [], [], opts);  AmBest = zeros(0, numel(var));
    assert(~isempty(c), 'No feasible parameter set (check boundFactor / limit / stability).');  return
end
costs = nan(2^nF, 1);  sols = cell(2^nF, 1);  rows = cell(2^nF, 1);
    function evalcode(code)
        if ~isnan(costs(code + 1)), return, end
        s = ones(nS, 1);  s(fq) = 1 - 2*bitget(code, 1:nF).';  Am = mono_rows(s, Nmin, nS, nc);  Am = Am(sort([2*fq(:); 2*fq(:) - 1]), var);   % free parameters, reduced columns
        cc = lsqlin(Rm, dm, [Ain; Am], [bin; zeros(size(Am, 1), 1)], Aeq, beq, lb, [], [], opts);
        if isempty(cc), costs(code + 1) = inf; else, costs(code + 1) = sum((Rm*cc - dm).^2);  sols{code + 1} = cc;  rows{code + 1} = Am; end
    end
if nF <= 6
    for code = 0:2^nF - 1, evalcode(code); end
else
    sp = ones(nS, 1);  sp(4:6) = -1;  starts = [0, 2^nF - 1, sum(2.^(find(sp(fq) < 0) - 1))];   % all decreasing, all increasing, couple increasing
    for code = starts
        evalcode(code);  improved = true;
        while improved
            improved = false;
            for b = 1:nF
                c2 = bitxor(code, 2^(b-1));  evalcode(c2);
                if costs(c2 + 1) < costs(code + 1), code = c2;  improved = true; end
            end
        end
    end
end
[best, ib] = min(costs);
assert(isfinite(best), 'No feasible parameter set (check boundFactor / limit / stability / monotone).');
c = sols{ib};  AmBest = rows{ib};
end

function Am = mono_rows(s, Nmin, nS, nc)
% monotone direction s_q (+1 decreasing, -1 increasing) for all N >= Nmin:  -s_q c1 <= 0  and  -s_q (c1 Nmin + 2 c2) <= 0
Am = zeros(2*nS, nc);
for q = 1:nS
    i1 = 3*(q-1) + 2;  Am(2*q-1, i1) = -s(q);  Am(2*q, [i1 i1+1]) = -s(q)*[Nmin, 2];
end
end

function [c, Jfinal, r] = refine(P, D, AJ, bJ, c, Ain, bin, Aeq, beq, lb, ub, var, cfix, nc, scale)
% Levenberg-Marquardt (damped Gauss-Newton) on  J = ||R1 c - d1||^2 + wd*||r||^2   (R1 c - d1: residual equations divided by 'scale';  r: forward errors), constraints as in the linear fit.
% c = reduced unknowns (free parameters and G); the full coefficient vector is  accumulate_c(c, var, cfix, nc)
mo = P.model;  Ns = P.ident.NList;  nN = numel(Ns);  nL = size(P.ident.trc, 3);
S = micromorphic_fe(P, 'ident', P.fe.nyCalib);                                       % coarser mesh, assembled once
dat.uf = zeros(nN, nL);  dat.ef = zeros(nN, nL);  dat.dom = false(nN, nL);           % full-scale dominant edge displacement and energy
for iN = 1:nN
    for l = 1:nL
        dat.dom(iN, l) = abs(D.het(iN).uy(l)) >= abs(D.het(iN).ux(l));
        if dat.dom(iN, l), dat.uf(iN, l) = D.het(iN).uy(l); else, dat.uf(iN, l) = D.het(iN).ux(l); end
        dat.ef(iN, l) = D.het(iN).twoU(l);
    end
end
[Qm, Rm] = qr(AJ, 0);  nb = scale;  R1 = Rm/nb;  d1 = (Qm.'*bJ)/nb;  wd = mo.dispWeight/(nN*nL);  n2 = nN*nL;
obj = @(c, r) sum((R1*c - d1).^2) + wd*sum(r.^2);
opts = optimoptions('lsqlin', 'Display', 'off');
fwd = @(cv) forward_reduced(cv, var, cfix, nc, S, P, dat);
[r, Jc] = fwd(c);
fprintf('   start:  J = %.4g   residual rows %.4f   max |error| displacement %.2f %%, energy %.2f %%\n', obj(c, r), norm(R1*c - d1), 100*max(abs(r(1:n2))), 100*max(abs(r(n2+1:end))));
lam = 1e-3;  n = numel(c);                                             % Levenberg-Marquardt damping (relative), adapted after every trial
for it = 1:mo.refineIter
    Ast = [R1; sqrt(wd)*Jc];  bst = [d1; sqrt(wd)*(Jc*c - r)];  J0 = obj(c, r);  ok = false;  mu0 = mean(diag(Ast.'*Ast));
    for trial = 1:12
        mu = lam*mu0;
        cn = lsqlin([Ast; sqrt(mu)*eye(n)], [bst; sqrt(mu)*c], Ain, bin, Aeq, beq, lb, ub, [], opts);       % min ||linearised J||^2 + mu ||c_new - c||^2, same constraints
        if ~isempty(cn) && all(isfinite(cn))
            rt = fwd(cn);
            if obj(cn, rt) < J0, ok = true;  break; end
        end
        lam = lam*10;
    end
    if ~ok, break; end
    c = cn;  lam = max(lam/5, 1e-9);  [r, Jc] = fwd(c);
    fprintf('   it %2d:  J = %.4g   residual rows %.4f   max |error| displacement %.2f %%, energy %.2f %%   (damping %.1e, G/C_hom = %.4g)\n', it, obj(c, r), norm(R1*c - d1), 100*max(abs(r(1:n2))), 100*max(abs(r(n2+1:end))), lam, c(end));
    if J0 - obj(c, r) < 1e-4*J0, break; end
end
Jfinal = obj(c, r);
end

function [r, J] = forward_reduced(cv, var, cfix, nc, S, P, dat)
% forward errors and Jacobian w.r.t. the reduced unknowns
c = accumulate_c(cv, var, cfix, nc);
if nargout > 1, [r, J] = forward_model(c, S, P, dat);  J = J(:, var); else, r = forward_model(c, S, P, dat); end
end

function c = accumulate_c(cv, var, cfix, nc)
% full coefficient vector from the reduced unknowns cv (free parameters and G) and the fixed values cfix
c = cfix;  c(var) = cv;
end

function [r, J] = forward_model(c, S, P, dat)
% relative errors  r = [(|u| - |u_full|)/|u_full| ; (2U - 2U_full)/2U_full]  (N-major, then load) of the dominant edge displacement and the energy of the
% forward micromorphic solution, and their Jacobian w.r.t. c (adjoint: du/dtheta = -lam' Kp U with K lam = e (edge functional), d(2U)/dtheta = -U' Kp U)
mo = P.model;  th0 = P.material.Chom(1, 1);  Ns = P.ident.NList;  nN = numel(Ns);  nL = size(P.ident.trc, 3);  Ly = P.ident.Ly;  nc = numel(c);  nd = size(S.F, 1);  nS = 9;
r = zeros(2*nN*nL, 1);  J = zeros(2*nN*nL, nc);  G = c(nc)*th0;  cred = reshape(c(1:3*nS), 3, nS).';
for iN = 1:nN
    N = Ns(iN);  phi = [1, 1/N, 1/N^2];  ell2 = (mo.ellFactor*Ly/N)^2;  th = th0*(cred*phi.');
    K = ell2*G*S.Kg;
    for p = 1:nS, K = K + th(p)*S.Kp{p}; end
    E = zeros(nd, nL);
    for l = 1:nL, if dat.dom(iN, l), E(:, l) = S.edgeY; else, E(:, l) = S.edgeX; end, end
    X = K(S.free, S.free) \ [S.F(S.free, :), E(S.free, :)];
    U = zeros(nd, nL);  U(S.free, :) = X(:, 1:nL);  Lam = zeros(nd, nL);  Lam(S.free, :) = X(:, nL+1:end);
    um = sum(E.*U, 1);  em = sum(S.F.*U, 1);  uf = dat.uf(iN, :);  ef = dat.ef(iN, :);  idx = (iN-1)*nL + (1:nL);
    r(idx) = ((abs(um) - abs(uf))./abs(uf)).';  r(nN*nL + idx) = ((em - ef)./ef).';
    if nargout > 1
        for p = 1:nS + 1
            if p <= nS, Kp = S.Kp{p};  sc = 1; else, Kp = S.Kg;  sc = ell2; end
            KU = Kp*U;  gu = (sign(um).*(-sc*sum(Lam.*KU, 1))./abs(uf)).';  ge = ((-sc*sum(U.*KU, 1))./ef).';
            if p <= nS
                for k = 1:3, J(idx, 3*(p-1) + k) = th0*phi(k)*gu;  J(nN*nL + idx, 3*(p-1) + k) = th0*phi(k)*ge; end
            else
                J(idx, nc) = th0*gu;  J(nN*nL + idx, nc) = th0*ge;
            end
        end
    end
end
end

function C = stiffness_at(c, N)
% [a b 0; b a 0; 0 0 d] at N from the coefficients c (3x3: rows a, b, d; columns c0, c1, c2)
th = c(:, 1) + c(:, 2)/N + c(:, 3)/N^2;
C = [th(1) th(2) 0; th(2) th(1) 0; 0 0 th(3)];
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

function Q = quad_setup(Id, rs)
% Gauss points (x, y) and weights W of the plate, bubble test functions (1 - s^2) P_i(s) and derivatives (reference coordinate s), edge points
[xg, wg] = gauss_legendre(rs.nQuad(1));  [yg, vg] = gauss_legendre(rs.nQuad(2));
Q.x = (xg + 1)/2*Id.Lx;  Q.y = (yg + 1)/2*Id.Ly;  Q.W = (vg(:)/2*Id.Ly)*(wg(:)/2*Id.Lx).';
[V, dV] = legendre_basis(xg, rs.test(1));  Q.FX = (1 - xg(:).^2).*V;  Q.dFX = -2*xg(:).*V + (1 - xg(:).^2).*dV;
[V, dV] = legendre_basis(yg, rs.test(2));  Q.FY = (1 - yg(:).^2).*V;  Q.dFY = -2*yg(:).*V + (1 - yg(:).^2).*dV;
[Q.xe, Q.we] = gauss_legendre(rs.nEdge);  Q.PXe = legendre_basis(Q.xe, rs.nTestEdge);
end

function [x, w] = gauss_legendre(n)
i = (1:n-1);  b = i./sqrt(4*i.^2 - 1);  [V, Dg] = eig(diag(b, 1) + diag(b, -1));  [x, k] = sort(diag(Dg));  w = 2*V(1, k).'.^2;
end

function F = fields_at(Ce, Cp, X, Y, Lx, Ly)
% data fields at the points (X, Y) (any shape; components in the 3rd dimension): E = eps_macro, P = psi, gam = E - P, derivatives of psi
n = numel(X);  [E, P, Px, Py] = deal(zeros(n, 3));
for k = 1:3
    E(:, k) = legendre_eval(Ce(:,:,k), X(:), Y(:), Lx, Ly);
    [P(:, k), Px(:, k), Py(:, k)] = legendre_eval(Cp(:,:,k), X(:), Y(:), Lx, Ly);
end
rs = @(A) reshape(A, [size(X) 3]);
F.E = rs(E);  F.P = rs(P);  F.gam = rs(E - P);  F.Px = rs(Px);  F.Py = rs(Py);
end

function [v, vx, vy] = legendre_eval(C, x, y, Lx, Ly)
% smooth field C ((degY+1) x (degX+1)) and its first derivatives at the points (x, y)
[Px, dPx] = legendre_basis(2*x(:)/Lx - 1, size(C, 2) - 1);  [Py, dPy] = legendre_basis(2*y(:)/Ly - 1, size(C, 1) - 1);
v = sum((Py*C).*Px, 2);  vx = sum((Py*C).*dPx, 2)*(2/Lx);  vy = sum((dPy*C).*Px, 2)*(2/Ly);
end

function Z = mulB(B, F)
% B*F at every point (components of F in the 3rd dimension)
Z = zeros(size(F));
for i = 1:3, for j = 1:3, if B(i, j) ~= 0, Z(:,:,i) = Z(:,:,i) + B(i, j)*F(:,:,j); end, end, end
end

function T = model_terms(F, ell2)
% terms belonging to every parameter p (1..3: a_g b_g d_g;  4..6: a_c b_c d_c;  7..9: a_p b_p d_p of C_psi;  10: G) with C = a B1 + b B2 + d B3,
% B1 = diag(1,1,0), B2 = [0 1 0;1 0 0;0 0 0], B3 = diag(0,0,1):
%   S{p}: stress     (C_gamma: B gam,  C_couple: B psi,  C_psi: 0)
%   M{p}: coefficient in dW/dpsi  (C_gamma: -B gam,  C_couple: B gam - B psi,  C_psi: B psi)
%   W{p}: energy density (C_gamma: 1/2 gam'B gam,  C_couple: psi'B gam,  C_psi: 1/2 psi'B psi),  W{10} = 1/2 l^2 |grad psi|^2 (times G)
B = {diag([1 1 0]), [0 1 0; 1 0 0; 0 0 0], diag([0 0 1])};
T.S = cell(1, 9);  T.M = cell(1, 9);  T.W = cell(1, 10);
for p = 1:3
    Bg = mulB(B{p}, F.gam);  Bp = mulB(B{p}, F.P);
    T.S{p}   = Bg;              T.M{p}   = -Bg;        T.W{p}   = 0.5*sum(F.gam.*Bg, 3);
    T.S{3+p} = Bp;              T.M{3+p} = Bg - Bp;    T.W{3+p} = sum(F.P.*Bg, 3);
    T.S{6+p} = zeros(size(Bp)); T.M{6+p} = Bp;         T.W{6+p} = 0.5*sum(F.P.*Bp, 3);
end
T.W{10} = 0.5*ell2*sum(F.Px.^2 + F.Py.^2, 3);
end

function blocks = residual_rows(Id, rs, mo, Chom, Ce, Cp, trc, Uhet, N, Q, xref)
% rows of  A*theta = b  (theta = [9 stiffness parameters, G]) for one N and one load, grouped by residual type (see the header)
Lx = Id.Lx;  Ly = Id.Ly;  th0 = Chom(1, 1);  ell2 = (mo.ellFactor*Ly/N)^2;  mx = rs.test(1);  my = rs.test(2);  nt = rs.nTestEdge;  w = rs.weights;  nP = 9;  nCol = nP + 1;
[X, Y] = meshgrid(Q.x, Q.y);  F = fields_at(Ce, Cp, X, Y, Lx, Ly);  T = model_terms(F, ell2);
blocks = struct('group', {}, 'A', {}, 'b', {});
% macro:  int sigma : eps(v) = 0,  v = (f,0) or (0,f)
A = zeros(0, nCol);
for c = 1:2
    for i = 1:mx+1
        for j = 1:my+1
            fx = (2/Lx)*Q.FY(:, j)*Q.dFX(:, i).';  fy = (2/Ly)*Q.dFY(:, j)*Q.FX(:, i).';
            if c == 1, e = {fx, 0*fx, fy}; else, e = {0*fx, fy, fx}; end
            row = zeros(1, nCol);
            for p = 1:nP, row(p) = sum(sum(Q.W.*(T.S{p}(:,:,1).*e{1} + T.S{p}(:,:,2).*e{2} + T.S{p}(:,:,3).*e{3}))); end
            A(end+1, :) = row; %#ok<AGROW>
        end
    end
end
blocks(end+1) = group('macro', A, zeros(size(A, 1), 1), xref, th0, w.macro);
% micro:  int (sum theta_p M_p) f + l^2 G grad(psi).grad(f) = 0
A = zeros(0, nCol);
for c = 1:3
    for i = 1:mx+1
        for j = 1:my+1
            f = Q.FY(:, j)*Q.FX(:, i).';  fx = (2/Lx)*Q.FY(:, j)*Q.dFX(:, i).';  fy = (2/Ly)*Q.dFY(:, j)*Q.FX(:, i).';
            row = zeros(1, nCol);
            for p = 1:nP, row(p) = sum(sum(Q.W.*T.M{p}(:,:,c).*f)); end
            row(nCol) = ell2*sum(sum(Q.W.*(F.Px(:,:,c).*fx + F.Py(:,:,c).*fy)));
            A(end+1, :) = row; %#ok<AGROW>
        end
    end
end
blocks(end+1) = group('micro', A, zeros(size(A, 1), 1), xref, th0, w.micro);
% traction on x = Lx (sigma_11 = tx, sigma_12 = ty) and free edges y = 0, Ly (sigma_12 = sigma_22 = 0)
ye = (Q.xe + 1)/2*Ly;  wye = Q.we/2*Ly;  xe = (Q.xe + 1)/2*Lx;  wxe = Q.we/2*Lx;
Te = model_terms(fields_at(Ce, Cp, Lx*ones(size(ye)), ye, Lx, Ly), ell2);
Pm = legendre_basis(2*ye/Ly - 1, 2);  tt = {Pm*trc(1, :).', [], Pm*trc(2, :).'};
A = zeros(0, nCol);  b = zeros(0, 1);
for k = 0:nt
    for comp = [1 3]
        row = zeros(1, nCol);
        for p = 1:nP, row(p) = sum(wye.*Te.S{p}(:,:,comp).*Q.PXe(:, k+1)); end
        A(end+1, :) = row;  b(end+1, 1) = sum(wye.*tt{comp}.*Q.PXe(:, k+1)); %#ok<AGROW>
    end
end
for yy = [0 Ly]
    Tf = model_terms(fields_at(Ce, Cp, xe, yy*ones(size(xe)), Lx, Ly), ell2);
    for k = 0:nt
        for comp = [3 2]
            row = zeros(1, nCol);
            for p = 1:nP, row(p) = sum(wxe.*Tf.S{p}(:,:,comp).*Q.PXe(:, k+1)); end
            A(end+1, :) = row;  b(end+1, 1) = 0; %#ok<AGROW>
        end
    end
end
blocks(end+1) = group('trac', A, b, xref, th0, w.trac);
% mu:  l^2 G (grad psi) n = 0 on the edges [x=Lx, x=0, y=Ly, y=0]
ex = {Lx*ones(size(ye)), 0*ye, xe, xe};  ey = {ye, ye, Ly*ones(size(xe)), 0*xe};  en = [1 0; -1 0; 0 1; 0 -1];  ew = {wye, wye, wxe, wxe};
A = zeros(0, nCol);
for q = find(rs.muEdges(:).' ~= 0)
    Fe = fields_at(Ce, Cp, ex{q}, ey{q}, Lx, Ly);
    for c = 1:3
        dn = en(q, 1)*Fe.Px(:,:,c) + en(q, 2)*Fe.Py(:,:,c);
        for k = 0:nt, row = zeros(1, nCol);  row(nCol) = ell2*sum(ew{q}.*dn.*Q.PXe(:, k+1));  A(end+1, :) = row; end %#ok<AGROW>
    end
end
blocks(end+1) = group('mu', A, zeros(size(A, 1), 1), xref, th0, w.mu);
% energy:  int W = U_het
row = zeros(1, nCol);
for p = 1:nCol, row(p) = sum(sum(Q.W.*T.W{p})); end
blocks(end+1) = group('energy', row, Uhet, xref, th0, w.energy);
% force and moment on the sections x = s Lx (equilibrium of the part x > s Lx)
Fx0 = Ly*trc(1, 1);  Fy0 = Ly*trc(2, 1);  Af = zeros(0, nCol);  bf = zeros(0, 1);  Am = zeros(0, nCol);  bm = zeros(0, 1);
for s = rs.sections
    xs = s*Lx;  Ts = model_terms(fields_at(Ce, Cp, xs*ones(size(ye)), ye, Lx, Ly), ell2);
    r1 = zeros(1, nCol);  r2 = r1;  r3 = r1;
    for p = 1:nP
        r1(p) = sum(wye.*Ts.S{p}(:,:,1));  r2(p) = sum(wye.*Ts.S{p}(:,:,3));  r3(p) = sum(wye.*(ye - Ly/2).*Ts.S{p}(:,:,1));
    end
    Af = [Af; r1; r2];  bf = [bf; Fx0; Fy0];  Am = [Am; r3];  bm = [bm; -((Lx - xs)*Fy0 - Ly^2*trc(1, 2)/6)]; %#ok<AGROW>
end
blocks(end+1) = group('force', Af, bf, xref, th0, w.force);
blocks(end+1) = group('moment', Am, bm, xref, th0, w.moment);
end

function g = group(name, A, b, xref, th0, w)
% unknowns in units of C_hom,11; rows scaled by the size of their terms; group weight w spread over its rows
if isempty(A), g = struct('group', name, 'A', A, 'b', b);  return; end
A = A*th0;  s = 1/(norm(abs(A)*xref) + norm(b) + eps);
g = struct('group', name, 'A', A*s*sqrt(w/size(A, 1)), 'b', b*s*sqrt(w/size(A, 1)));
end
