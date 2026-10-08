function C = compare_results(P, D, M)
% =========================================================================
%  COMPARISON  micromorphic theory (M)  vs  full-scale plate (D, centred inclusion)  vs  homogeneous plate         C = compare_results(P, D, M)
% =========================================================================
%  For every N and load (errors in %, > 0: more flexible than the full-scale plate):
%      errDisp   (|u_micro| - |u_full|)/|u_full|     dominant mean displacement of the edge x = Lx            (errHomDisp: homogeneous plate)
%      errEn     (2U_micro - 2U_full)/2U_full         2U = F'U                                                  (errHomEn)
%      errPsi    ||psi_micro - eps_full||/||eps_full||   inside the plate (15% margins removed)
%      errSig    ||sigma_micro - sigma_full||/||sigma_full||                                                     (errSigHom: homogeneous plate)
%      secF, secM  rms error of the section force  (int sigma_11 dy, int sigma_12 dy)  and moment  (int (y-Ly/2) sigma_11 dy)  of the stress
%                  against the applied load, on interior sections x = const                          (.micro, .full, .hom)
%      errReact  clamp reaction R = K U - F against the exact value -(applied load)  (check of the FE implementation)
%  All results are [nN x nLoads] arrays.
% =========================================================================
T = P.test;  Lx = T.Lx;  Ly = T.Ly;  Ns = T.NList;  nN = numel(Ns);  nL = size(T.trc, 3);  names = T.loads;  THK = P.material.thickness;
xc = M.xc;  yc = M.yc;  dB = M.dB;
mask  = (((yc/Ly > 0.15) & (yc/Ly < 0.85)).') * ((xc/Lx > 0.15) & (xc/Lx < 0.95)) > 0;      % interior bins for the field misfits
xmask = (xc/Lx > 0.10) & (xc/Lx < 0.90);                                                      % interior sections
z = zeros(nN, nL);  C = struct('Ns', Ns, 'names', {names}, 'refDom', z, 'mmDom', z, 'homDom', z, 'errDisp', z, 'errHomDisp', z, 'errEn', z, 'errHomEn', z, ...
    'errPsi', z, 'errSig', z, 'errSigHom', z, 'secF', struct('micro', z, 'full', z, 'hom', z), 'secM', struct('micro', z, 'full', z, 'hom', z), ...
    'errReact', struct('micro', z, 'full', z, 'hom', z));
for iN = 1:nN
    f = D.het(iN);  m = M.N(iN);
    for l = 1:nL
        if abs(f.uy(l)) >= abs(f.ux(l)), rv = f.uy(l);  mv = m.uy(l);  hv = D.hom.uy(l); else, rv = f.ux(l);  mv = m.ux(l);  hv = D.hom.ux(l); end
        C.refDom(iN, l) = rv;  C.mmDom(iN, l) = mv;  C.homDom(iN, l) = hv;
        C.errDisp(iN, l) = (abs(mv) - abs(rv))/abs(rv)*100;  C.errHomDisp(iN, l) = (abs(hv) - abs(rv))/abs(rv)*100;
        C.errEn(iN, l) = (m.twoU(l) - f.twoU(l))/f.twoU(l)*100;  C.errHomEn(iN, l) = (D.hom.twoU(l) - f.twoU(l))/f.twoU(l)*100;
        Pf = grid_eval(f.psiCoef(:,:,:,l), xc, yc, Lx, Ly);  Sf = grid_eval(f.sigCoef(:,:,:,l), xc, yc, Lx, Ly);  Sh = grid_eval(D.hom.sigCoef(:,:,:,l), xc, yc, Lx, Ly);
        Sm = m.Sbin(:,:,:,l);
        C.errPsi(iN, l) = misfit(m.Pbin(:,:,:,l), Pf, mask);  C.errSig(iN, l) = misfit(Sm, Sf, mask);  C.errSigHom(iN, l) = misfit(Sh, Sf, mask);
        tl = T.trc(:,:,l);  Fx0 = Ly*tl(1, 1);  Fy0 = Ly*tl(2, 1);  M0 = -((Lx - xc)*Fy0 - Ly^2*tl(1, 2)/6);       % applied load resultants
        Fsc = Ly*max(max(abs(tl(:))), 1e-12);  Msc = max(max(abs(M0)), Ly^2*max(abs(tl(:)))/6);
        S = {Sm, Sf, Sh};  fn = {'micro', 'full', 'hom'};
        for k = 1:3
            R = [sum(S{k}(:,:,1), 1)*dB; sum(S{k}(:,:,3), 1)*dB; sum(S{k}(:,:,1).*repmat(yc(:) - Ly/2, 1, numel(xc)), 1)*dB];
            C.secF.(fn{k})(iN, l) = 100*sqrt(mean((R(1, xmask) - Fx0).^2 + (R(2, xmask) - Fy0).^2))/Fsc;
            C.secM.(fn{k})(iN, l) = 100*sqrt(mean((R(3, xmask) - M0(xmask)).^2))/Msc;
        end
        ex = -THK*[Ly*tl(1, 1); Ly*tl(2, 1); Lx*Ly*tl(2, 1) - Ly^2*tl(1, 2)/6];  sc = max(max(abs(ex)), eps);
        C.errReact.micro(iN, l) = 100*max(abs(m.react(:, l) - ex))/sc;  C.errReact.full(iN, l) = 100*max(abs(f.react(:, l) - ex))/sc;
        C.errReact.hom(iN, l) = 100*max(abs(D.hom.react(:, l) - ex))/sc;
    end
end
for l = 1:nL
    fprintf('\n  load ''%s'': errors in %% w.r.t. the full-scale plate (> 0: more flexible)\n', names{l});
    fprintf('  %3s %7s | %12s %12s %12s | %8s %8s | %8s %8s | %7s %8s %8s | %8s %8s\n', 'N', 'cell', 'full-scale', 'homogeneous', 'micromorph.', 'd micro', 'd hom', '2U micro', '2U hom', 'psi', 'sig mic', 'sig hom', 'F micro', 'M micro');
    for iN = 1:nN
        fprintf('  %3d %7.4f | %12.5e %12.5e %12.5e | %+7.2f%% %+7.2f%% | %+7.2f%% %+7.2f%% | %6.1f%% %7.1f%% %7.1f%% | %7.2f%% %7.2f%%\n', Ns(iN), D.het(iN).Lcell, ...
            C.refDom(iN, l), C.homDom(iN, l), C.mmDom(iN, l), C.errDisp(iN, l), C.errHomDisp(iN, l), C.errEn(iN, l), C.errHomEn(iN, l), ...
            C.errPsi(iN, l), C.errSig(iN, l), C.errSigHom(iN, l), C.secF.micro(iN, l), C.secM.micro(iN, l));
    end
    fprintf('  FE clamp reaction vs exact (%%): micro %s, full %s, hom %s\n', mat2str(C.errReact.micro(:, l).', 2), mat2str(C.errReact.full(:, l).', 2), mat2str(C.errReact.hom(:, l).', 2));
end
mx = max(abs(C.errDisp(:)));
fprintf('\n  largest |displacement error| of the micromorphic solution: %.2f %% -> %s (target %.1f %%);  homogeneous: %.2f %%\n', mx, ternary(mx <= P.out.tolerance, 'PASS', 'NOT reached'), P.out.tolerance, max(abs(C.errHomDisp(:))));
fprintf('  largest |energy error| of the micromorphic solution: %.2f %%;  homogeneous: %.2f %%\n', max(abs(C.errEn(:))), max(abs(C.errHomEn(:))));
end

function e = misfit(A, B, mask)
A3 = reshape(A, [], 3);  B3 = reshape(B, [], 3);  mk = mask(:);
e = 100*norm(A3(mk, :) - B3(mk, :), 'fro')/norm(B3(mk, :), 'fro');
end

function Z = grid_eval(C, xc, yc, Lx, Ly)
% smooth Legendre fields C ((degY+1) x (degX+1) x 3) on the tensor grid of the bin centres -> numel(yc) x numel(xc) x 3
Px = legendre_basis(2*xc(:)/Lx - 1, size(C, 2) - 1);  Py = legendre_basis(2*yc(:)/Ly - 1, size(C, 1) - 1);
Z = zeros(numel(yc), numel(xc), 3);
for k = 1:3, Z(:, :, k) = Py*C(:, :, k)*Px.'; end
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end
