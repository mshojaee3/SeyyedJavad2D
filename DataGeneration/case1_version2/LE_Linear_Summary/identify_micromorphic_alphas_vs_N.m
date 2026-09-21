% ============================================================
%  identify_micromorphic_alphas_vs_N.m
%
%  Reads every Ensemble_Summary_NxN.csv (+ its _homogenized.csv pair)
%  found in THIS SCRIPT'S OWN FOLDER, identifies the micromorphic
%  parameters under the ansatz
%       sym(psi)(X) = H + sym(G.(X-Xc))          [Eq. psi=H+G.X]
%  with alpha_m = 1 FIXED (Affine Equilibrium Ansatz / Cauchy-
%  consistency requirement -- see project notes), and plots the
%  identified parameters vs N, plus fit-quality diagnostics.
%
%  NOTATION: everything below is written directly in terms of
%  sym(psi)(X) and eps_M (never introducing the relative strain
%  gamma=eps_M-sym(psi) as a named quantity), since sym(psi) is
%  position-dependent while eps_M is a single value -- gamma would
%  only ever appear as their pointwise difference, which is what J
%  below already is.
%
%  WHAT IS AND ISN'T IDENTIFIABLE FROM THIS DATA -- READ THIS FIRST
%  ------------------------------------------------------------
%  int_V sym(psi)(X) dV = eps_M*V exactly (the G-part of the ansatz
%  integrates to zero over a domain symmetric about its own center).
%  Substituting sym(psi)=eps_M-... directly into W and integrating,
%  using J := int_V (sym(psi)-eps_M)'*C*(sym(psi)-eps_M) dV (the
%  variance of sym(psi) about its own domain mean -- computable in
%  closed form from G, C and the domain size, no gamma needed), gives
%       int_V W dV = (1/2)*V*eps_M'*C*eps_M + (1+alpha_h-2*alpha_c)/2 * J.
%  So the domain-INTEGRATED energy depends on alpha_h and alpha_c ONLY
%  through the single combination
%       beta := 1 + alpha_h - 2*alpha_c
%  for ANY load case -- adding more load cases does NOT break this,
%  it is structural to this ansatz (matches the "rank deficiency"
%  already noted in the 1D-bar project), not a data limitation.
%
%  alpha_ell IS separately identifiable: the double stress
%  mu = d W/d kappa = alpha_ell*(L/N)^2*A_kappa*kappa depends on
%  alpha_ell ONLY (not alpha_h, alpha_c), so Q-matching gives alpha_ell
%  cleanly.
%
%  So this script reports, per N:
%    - beta        (well-posed, from ENERGY:  ALLSE_het - ALLSE_macro)
%    - alpha_ell   (well-posed, from Q:        Q_het - Q_hom)
%    - alpha_h, alpha_c individually -- via the MINIMUM-NORM point on
%      the degenerate line alpha_h-2*alpha_c=beta-1. This is a
%      reproducible NUMBER, not a genuine separate identification --
%      it is plotted only because it was asked for; treat it as
%      illustrative until psi is identified independently of (H,G)
%      (e.g. via the windowed least-squares field projection used in
%      the 1D-bar project), which is the only way to break the
%      degeneracy with this kind of data.
%
%  PLUS, for every N: parity plots (true vs model, expect y=x) and
%  RMSE / R^2 / residual-histogram diagnostics for both the energy
%  fit and the Q fit, using the SAME beta, alpha_ell just identified
%  for that N.
%
%  ENERGY REGRESSOR: J_meas, not the ansatz-theoretical J -- IMPORTANT
%  -------------------------------------------------------------------
%  The theoretical J = int_V(sym(psi)-eps_M)'C(sym(psi)-eps_M)dV (call
%  it J_theory) assumes psi=H+G.X is a good approximation of the TRUE
%  local strain field in a homogeneous medium under this same G. In
%  practice this assumption is bad: the diagnostic below routinely
%  shows |ALLSE_hom-(macro+J_theory)| of the SAME ORDER as ALLSE
%  itself, not a small correction. Since J_theory is what the energy
%  fit regresses beta against, a badly-wrong J_theory shows up
%  directly as scatter in the energy parity plot -- this is what you
%  are likely seeing if that plot looks noisy rather than tight on
%  y=x. The fix used below: replace J_theory with the MEASURED
%       J_meas(case) := ALLSE_hom(case) - ALLSE_macro(case)
%  i.e. the actual FE homogeneous-solve energy for that case's own G,
%  which needs no assumption about psi matching the true field at all.
%  J_theory is still computed and reported (diag_energy_ansatz =
%  J_meas - J_theory) purely as a diagnostic of how bad that ansatz
%  assumption is -- it no longer enters the fit.
%
%  BASELINE VALUES: beta -> 2, not 0, as N -> infinity -- IMPORTANT
%  -------------------------------------------------------------------
%  J_meas is a FIXED, N-independent classical-elasticity quantity (no
%  microstructure involved). Since E_micro = ALLSE_het-ALLSE_macro ->
%  J_meas as N->infinity (the KUBC finite-size limit: ALLSE_het ->
%  ALLSE_hom), and E_micro = (beta/2)*J_meas + alpha_ell-term, beta
%  itself must tend to 2 (not 0) as genuine microstructure vanishes.
%  betaMinus2 := beta-2 is the quantity that actually goes to zero --
%  it is reported and plotted alongside raw beta.
%  Similarly, raw alpha_ell need NOT go to zero even as the true
%  microstructural signal vanishes: it is defined via an ASSUMED
%  ellc^2 = (L/N)^2 power law, and if the true finite-size decay in
%  Q_het-Q_hom is closer to 1/N (as it empirically is for energy),
%  alpha_ell = (that 1/N signal)/ellc^2 will grow like N to
%  compensate. muScale := alpha_ell*ellc^2 (the raw fitted slope of
%  Q_micro vs kappa, with no assumed power of ellc) is the more
%  reliable quantity to watch for convergence.
%
%  STRESS: the macro/Cauchy stress is NOT used in the objective --
%  per the project's own established finding, it depends only on
%  alpha_m (fixed=1) and C_hom (already known), so it carries zero
%  information about alpha_h, alpha_c, alpha_ell once those are fixed.
%  It is reported as a pure consistency check instead (C_hom*H vs
%  the homogenized file's own Savg -- should agree to ~machine
%  precision; a large mismatch would flag a C_hom or data problem).
%
%  DEPENDENCY: homogenize2D_PBC.m must be on the MATLAB path or in
%  this same folder (delivered earlier in this pipeline).
%
%  NOTE: written carefully but NOT executed here (no MATLAB in this
%  environment) -- please run it and report back any error. The two
%  places most worth double-checking against your own conventions are
%  (a) the A_kappa weighting (defaulted to identity below -- replace
%  if your paper fixes specific weights) and (b) the Q111..Q222 ->
%  6-component kappa mapping (derived from your own domainAverages_LE
%  code and printed as a sanity check below).
% ============================================================

function identify_micromorphic_alphas_vs_N()

scriptDir = fileparts(mfilename('fullpath'));

%% ============ USER INPUT ============
% Material/geometry -- MUST match your main RUC pipeline exactly: these
% determine C_hom via homogenize2D_PBC.m, independent of N.
mat.E_m = 70000;  mat.nu_m = 0.33;
mat.E_i = 3500;   mat.nu_i = 0.33;
geom.Rfrac = 6/19; geom.meshFrac = 0.02; geom.N = 1;
geom.Lx_tot = 1.0; geom.Ly_tot = 1.0;

Lx_tot = geom.Lx_tot; Ly_tot = geom.Ly_tot;

% A_kappa = diag(a1..a6) for kappa = [kxxx,kxxy,kyyx,kyyy,kxyx,kxyy].
% NOT specified anywhere in the CSV data -- defaults to identity.
% Replace with your paper's actual weights if they differ.
Akappa_diag = ones(1,6);

OUTNAME = fullfile(scriptDir, 'Alphas_vs_N.csv');
DETAILNAME = fullfile(scriptDir, 'Alphas_vs_N_detail.csv');

%% ============ DISCOVER N's PRESENT IN THIS FOLDER ============
files = dir(fullfile(scriptDir, 'Ensemble_Summary_*.csv'));
Ns = [];
for k = 1:numel(files)
    fn = files(k).name;
    if contains(fn, 'homogenized'), continue; end
    tok = regexp(fn, '^Ensemble_Summary_(\d+)X(\d+)\.csv$', 'tokens', 'once');
    if isempty(tok), continue; end
    n1 = str2double(tok{1}); n2 = str2double(tok{2});
    if n1 ~= n2
        warning('Skipping %s: NX (%d) ~= NY (%d), this script assumes square NxN.', fn, n1, n2);
        continue;
    end
    Ns(end+1) = n1; %#ok<AGROW>
end
Ns = sort(unique(Ns));
if isempty(Ns)
    error('identify_micromorphic_alphas_vs_N:NoFiles', ...
        'No Ensemble_Summary_NxN.csv files found in %s.', scriptDir);
end
fprintf('Found N = %s in %s\n\n', mat2str(Ns), scriptDir);

%% ============ C_hom, computed ONCE (independent of N) ============
if exist('homogenize2D_PBC', 'file') ~= 2
    error('identify_micromorphic_alphas_vs_N:MissingDep', ...
        'homogenize2D_PBC.m not found on the MATLAB path or in this folder.');
end
[C_hom, Cs33_hom] = homogenize2D_PBC(mat, geom); %#ok<ASGLU>
fprintf('C_hom (Voigt, [S11;S22;S12]=C*[E11;E22;gamma12]) =\n'); disp(C_hom);

%% ============ PER-N IDENTIFICATION ============
nN = numel(Ns);
results = struct('N',{},'ellc',{},'beta',{},'betaMinus2',{},'alphaL',{},'muScale',{}, ...
                  'alphaH_minnorm',{},'alphaC_minnorm',{}, ...
                  'diag_energy_ansatz',{},'diag_stress_consistency',{}, ...
                  'RMSE_E',{},'RMSE_Q',{},'R2_E',{},'R2_Q',{}, ...
                  'nCases',{});
detN = []; detCase = {}; detEtrue = []; detEmodel = [];
detQtrue = zeros(0,6); detQmodel = zeros(0,6);

for ii = 1:nN
    N = Ns(ii);
    fHet = fullfile(scriptDir, sprintf('Ensemble_Summary_%dX%d.csv', N, N));
    fHom = fullfile(scriptDir, sprintf('Ensemble_Summary_%dX%d_homogenized.csv', N, N));
    if ~isfile(fHom)
        warning('Skipping N=%d: %s not found.', N, fHom);
        continue;
    end
    Thet = readtable(fHet);
    Thom = readtable(fHom);

    [commonCases, ia, ib] = intersect(Thet.LoadCase, Thom.LoadCase, 'stable');
    nCases = numel(commonCases);
    if nCases == 0
        warning('Skipping N=%d: no matching LoadCase between the two files.', N);
        continue;
    end

    ellc = Lx_tot / N;   % unit-cell length for this N (ell_c = L/N)

    Iee_all       = zeros(nCases,1);   % J_theory (diagnostic only, see header)
    Jmeas_all     = zeros(nCases,1);   % J_meas = ALLSE_hom(case)-ALLSE_macro(case): the energy-fit regressor
    kappaAk2_all  = zeros(nCases,1);   % kappa' * (Akappa .* kappa), per case
    Emicro_all    = zeros(nCases,1);   % ALLSE_het - ALLSE_macro(closed form)
    V_all         = zeros(nCases,1);
    diagE_all     = zeros(nCases,1);   % ALLSE_hom - (ALLSE_macro + J) [diagnostic]
    stressDiag_all= zeros(nCases,1);   % max|C_hom*H - Savg_hom|         [diagnostic]
    kappa_all     = zeros(nCases,6);   % stored per-case for later model-prediction plots
    Qmicro_all    = zeros(nCases,6);   % Q_het - Q_hom, per case
    regQ_v = []; regQ_t = [];

    for c = 1:nCases
        rh = Thet(ia(c),:);
        rH = Thom(ib(c),:);

        % ---- H, from the (machine-precision-recovered) Fbar ----
        H  = [rh.Fbar_11-1, rh.Fbar_12; rh.Fbar_12, rh.Fbar_22-1];
        Hv = [H(1,1); H(2,2); 2*H(1,2)];

        % ---- kappa (6 indep. components) from G111..G222 (Method A) ----
        % mapping: kappa = [k_xxx,k_xxy,k_yyx,k_yyy,k_xyx,k_xyy]
        %                 = [G111, G112, G221, G222, 0.5(G121+G211), 0.5(G122+G212)]
        kappa = [rh.G111; rh.G112; rh.G221; rh.G222; ...
                 0.5*(rh.G121+rh.G211); 0.5*(rh.G122+rh.G212)];
        Pv = [kappa(1); kappa(3); 2*kappa(5)];   % Voigt coeff. of (X1-Xc1)
        Qv = [kappa(2); kappa(4); 2*kappa(6)];   % Voigt coeff. of (X2-Xc2)

        V = rh.Volume;
        Iee = V*((Lx_tot^2/12)*(Pv.'*C_hom*Pv) + (Ly_tot^2/12)*(Qv.'*C_hom*Qv));
        ALLSE_macro = 0.5*V*(Hv.'*C_hom*Hv);
        Jmeas = rH.ALLSE - ALLSE_macro;   % MEASURED curvature energy for this case's own G

        Iee_all(c)      = Iee;
        Jmeas_all(c)    = Jmeas;
        V_all(c)        = V;
        kappaAk2_all(c) = kappa.' * (Akappa_diag(:).*kappa);
        Emicro_all(c)   = rh.ALLSE - ALLSE_macro;

        % ---- diagnostic (NOT used in the fit): how badly does the
        %      psi=H+G.X ansatz's own theoretical curvature energy
        %      (Iee) predict the ACTUAL homogeneous energy excess
        %      (Jmeas)? A large value here is exactly why the fit
        %      below uses Jmeas, not Iee. ----
        diagE_all(c) = Jmeas - Iee;
        Sv_hom = [rH.Savg_11; rH.Savg_22; rH.Savg_12];
        stressDiag_all(c) = max(abs(C_hom*Hv - Sv_hom));

        % ---- Q target: heterogeneous excess over the homogeneous
        %      (purely-geometric, non-microstructural) baseline ----
        Qhet6 = [rh.Q111; rh.Q112; rh.Q221; rh.Q222; rh.Q121; rh.Q122];
        Qhom6 = [rH.Q111; rH.Q112; rH.Q221; rH.Q222; rH.Q121; rH.Q122];
        Qmicro6 = Qhet6 - Qhom6;

        kappa_all(c,:)  = kappa.';
        Qmicro_all(c,:) = Qmicro6.';

        regQ_v = [regQ_v; ellc^2*(Akappa_diag(:).*kappa)]; %#ok<AGROW>
        regQ_t = [regQ_t; Qmicro6]; %#ok<AGROW>
    end

    % ---- alpha_ell: closed-form scalar least squares over all
    %      (case, kappa-component) pairs stacked together ----
    if norm(regQ_v) > 0
        alphaL = (regQ_v.'*regQ_t) / (regQ_v.'*regQ_v);
    else
        alphaL = NaN;
    end

    % ---- beta = 1+alpha_h-2*alpha_c: closed-form scalar least
    %      squares over cases, given alpha_ell from the Q-fit --
    %      regressed against J_meas (measured), NOT the ansatz-
    %      theoretical Iee -- see header note ----
    rhs = Emicro_all - alphaL*(ellc^2/2).*kappaAk2_all;
    if norm(Jmeas_all) > 0
        beta = 2*sum(rhs.*Jmeas_all) / sum(Jmeas_all.^2);
    else
        beta = NaN;
    end

    % ---- betaMinus2: the quantity that ACTUALLY vanishes as N->infty
    %      (see header note) -- beta itself tends to 2, not 0, because
    %      J_meas is a fixed, N-independent classical-elasticity term.
    betaMinus2 = beta - 2;

    % ---- muScale: the raw fitted slope of Q_micro vs kappa, with NO
    %      assumption about which power of ellc it should scale with --
    %      this is what should shrink toward zero as N grows, even
    %      though alphaL itself (muScale/ellc^2) need not. ----
    muScale = alphaL*ellc^2;

    % ---- minimum-norm split of the degenerate line
    %      alpha_h - 2*alpha_c = beta-1  (see header note: NOT a
    %      genuine separate identification, just a reproducible point) ----
    k_ = beta - 1;
    alphaH_mn = 0.2*k_;
    alphaC_mn = -0.4*k_;

    % ---- per-case model predictions, using THIS N's just-identified
    %      beta, alphaL -- for the parity plots / RMSE / R^2 below ----
    Emodel_all = (beta/2)*Jmeas_all + alphaL*(ellc^2/2)*V_all.*kappaAk2_all;
    Qmodel_all = alphaL*ellc^2*(kappa_all .* Akappa_diag(:).');   % nCases x 6

    Eresid = Emicro_all - Emodel_all;
    Qresid = Qmicro_all - Qmodel_all;

    RMSE_E = sqrt(mean(Eresid.^2));
    RMSE_Q = sqrt(mean(Qresid(:).^2));
    SStotE = sum((Emicro_all - mean(Emicro_all)).^2);
    SStotQ = sum((Qmicro_all(:) - mean(Qmicro_all(:))).^2);
    R2_E = 1 - sum(Eresid.^2)/max(SStotE, eps);
    R2_Q = 1 - sum(Qresid(:).^2)/max(SStotQ, eps);

    detN      = [detN; repmat(N, nCases, 1)]; %#ok<AGROW>
    detCase   = [detCase; commonCases(:)]; %#ok<AGROW>
    detEtrue  = [detEtrue; Emicro_all]; %#ok<AGROW>
    detEmodel = [detEmodel; Emodel_all]; %#ok<AGROW>
    detQtrue  = [detQtrue; Qmicro_all]; %#ok<AGROW>
    detQmodel = [detQmodel; Qmodel_all]; %#ok<AGROW>

    results(end+1) = struct('N',N, 'ellc',ellc, 'beta',beta, 'betaMinus2',betaMinus2, ...
        'alphaL',alphaL, 'muScale',muScale, ...
        'alphaH_minnorm',alphaH_mn, 'alphaC_minnorm',alphaC_mn, ...
        'diag_energy_ansatz',mean(diagE_all), ...
        'diag_stress_consistency',max(stressDiag_all), ...
        'RMSE_E',RMSE_E, 'RMSE_Q',RMSE_Q, 'R2_E',R2_E, 'R2_Q',R2_Q, ...
        'nCases',nCases); %#ok<AGROW>

    fprintf(['N=%2d  ellc=%.4f  nCases=%d | beta=1+aH-2aC=%8.4f  alphaL=%8.4f | ' ...
             '(min-norm) alphaH=%7.4f alphaC=%7.4f | RMSE_E=%.3e R2_E=%.4f  RMSE_Q=%.3e R2_Q=%.4f | ' ...
             'diag: Jmeas-Jtheory=%.3e (not used in fit)  |C*H-S|=%.3e\n'], ...
        N, ellc, nCases, beta, alphaL, alphaH_mn, alphaC_mn, RMSE_E, R2_E, RMSE_Q, R2_Q, ...
        results(end).diag_energy_ansatz, results(end).diag_stress_consistency);
end

if isempty(results)
    error('identify_micromorphic_alphas_vs_N:NoValidN', 'No N had a valid matching file pair.');
end

%% ============ WRITE SUMMARY + DETAIL CSVs ============
Tout = struct2table(results);
writetable(Tout, OUTNAME);
fprintf('\nWrote %s\n', OUTNAME);

Qnames_true  = {'Q111_true','Q112_true','Q221_true','Q222_true','Q121_true','Q122_true'};
Qnames_model = {'Q111_model','Q112_model','Q221_model','Q222_model','Q121_model','Q122_model'};
Tdetail = table(detN, detCase, detEtrue, detEmodel, 'VariableNames', {'N','LoadCase','E_true','E_model'});
Tdetail = [Tdetail, array2table(detQtrue,'VariableNames',Qnames_true), ...
                     array2table(detQmodel,'VariableNames',Qnames_model)];
writetable(Tdetail, DETAILNAME);
fprintf('Wrote %s\n', DETAILNAME);

%% ============ PLOTS ============
Nv = Tout.N;

figure('Color','w','Name','beta and alpha_ell vs N');
subplot(1,2,1);
plot(Nv, Tout.beta, '-o', 'LineWidth',1.5, 'MarkerSize',6);
grid on; xlabel('N'); ylabel('\beta = 1+\alpha_h-2\alpha_c');
title('Well-posed: \beta vs N (from energy)');

subplot(1,2,2);
plot(Nv, Tout.alphaL, '-s', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0.85 0.33 0.10]);
grid on; xlabel('N'); ylabel('\alpha_\ell');
title('Well-posed: \alpha_\ell vs N (from Q)');

figure('Color','w','Name','beta-2 and mu_scale: the quantities that should -> 0 at large N');
subplot(1,2,1);
plot(Nv, Tout.betaMinus2, '-o', 'LineWidth',1.5, 'MarkerSize',6);
grid on; xlabel('N'); ylabel('\beta - 2');
title({'\beta \rightarrow 2 (not 0) as N\rightarrow\infty is EXPECTED --', ...
       'J_{meas} is a fixed classical term. \beta-2 is what vanishes.'}, ...
       'FontWeight','normal','FontSize',9);

subplot(1,2,2);
plot(Nv, Tout.muScale, '-s', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0.85 0.33 0.10]);
grid on; xlabel('N'); ylabel('\mu_{scale} = \alpha_\ell \cdot \ell_c^2');
title({'Raw fitted Q_{micro}-vs-\kappa slope, no assumed \ell_c power --', ...
       'watch THIS for convergence, not raw \alpha_\ell (see header note)'}, ...
       'FontWeight','normal','FontSize',9);

figure('Color','w','Name','alpha_h, alpha_c (minimum-norm, illustrative) vs N');
plot(Nv, Tout.alphaH_minnorm, '-o', Nv, Tout.alphaC_minnorm, '-s', 'LineWidth',1.5, 'MarkerSize',6);
grid on; xlabel('N'); ylabel('\alpha');
legend({'\alpha_h (min-norm)','\alpha_c (min-norm)'}, 'Location','best');
title({'NOT a genuine separate identification -- see header note.', ...
       'Minimum-norm point on the degenerate line \alpha_h-2\alpha_c=\beta-1.'}, ...
       'FontWeight','normal', 'FontSize',9);

figure('Color','w','Name','Diagnostics vs N');
subplot(1,2,1);
semilogy(Nv, abs(Tout.diag_energy_ansatz), '-o', 'LineWidth',1.5, 'MarkerSize',6);
grid on; xlabel('N'); ylabel('|J_{meas} - J_{theory}|');
title({'psi-ansatz check (NOT used in the fit -- the fit uses J_{meas})', ...
       'Large values mean the psi=H+G.X field is a poor stand-in for the true one'}, ...
       'FontWeight','normal', 'FontSize',9);

subplot(1,2,2);
semilogy(Nv, Tout.diag_stress_consistency, '-o', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0.47 0.67 0.19]);
grid on; xlabel('N'); ylabel('max|C_{hom}H - S_{avg,hom}|');
title('C_{hom} / KUBC consistency check (should be ~machine eps)');

%% ---- Parity plot: energy, true vs model, colored by N, expect y=x ----
figure('Color','w','Name','Parity plot: energy');
cmap = lines(numel(Nv));
hold on;
hS = gobjects(numel(Nv),1);
for ii = 1:numel(Nv)
    m = (detN == Nv(ii));
    hS(ii) = scatter(detEtrue(m), detEmodel(m), 18, cmap(ii,:), 'filled', 'MarkerFaceAlpha',0.6);
end
lims = [min([detEtrue;detEmodel]), max([detEtrue;detEmodel])];
plot(lims, lims, 'k--', 'LineWidth',1.2);
hold off; grid on; axis equal; xlim(lims); ylim(lims);
xlabel('E_{micro} true  (ALLSE_{het} - ALLSE_{macro})');
ylabel('E_{micro} model  (\beta/2 \cdot J_{meas} + \alpha_\ell term)');
title('Energy: true vs model (dashed = y=x)');
legend(hS, "N="+string(Nv), 'Location','bestoutside');

%% ---- Parity plot: double stress Q, all 6 components, colored by N ----
figure('Color','w','Name','Parity plot: double stress Q');
hold on;
hS2 = gobjects(numel(Nv),1);
for ii = 1:numel(Nv)
    m = (detN == Nv(ii));
    hS2(ii) = scatter(reshape(detQtrue(m,:),[],1), reshape(detQmodel(m,:),[],1), 14, cmap(ii,:), 'filled', 'MarkerFaceAlpha',0.5);
end
limsQ = [min([detQtrue(:);detQmodel(:)]), max([detQtrue(:);detQmodel(:)])];
plot(limsQ, limsQ, 'k--', 'LineWidth',1.2);
hold off; grid on; axis equal; xlim(limsQ); ylim(limsQ);
xlabel('Q_{micro} true  (Q_{het}-Q_{hom}), all 6 components');
ylabel('Q_{micro} model  (\alpha_\ell \ell_c^2 A_\kappa \kappa)');
title('Double stress: true vs model (dashed = y=x)');
legend(hS2, "N="+string(Nv), 'Location','bestoutside');

%% ---- RMSE and R^2 vs N, for both fits ----
figure('Color','w','Name','Fit quality vs N');
subplot(2,2,1);
semilogy(Nv, Tout.RMSE_E, '-o', 'LineWidth',1.5, 'MarkerSize',6);
grid on; xlabel('N'); ylabel('RMSE'); title('Energy fit: RMSE vs N');
subplot(2,2,2);
plot(Nv, Tout.R2_E, '-o', 'LineWidth',1.5, 'MarkerSize',6);
grid on; xlabel('N'); ylabel('R^2'); ylim([min(0,min(Tout.R2_E))-0.05, 1.05]);
title('Energy fit: R^2 vs N');
subplot(2,2,3);
semilogy(Nv, Tout.RMSE_Q, '-s', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0.85 0.33 0.10]);
grid on; xlabel('N'); ylabel('RMSE'); title('Q fit: RMSE vs N');
subplot(2,2,4);
plot(Nv, Tout.R2_Q, '-s', 'LineWidth',1.5, 'MarkerSize',6, 'Color',[0.85 0.33 0.10]);
grid on; xlabel('N'); ylabel('R^2'); ylim([min(0,min(Tout.R2_Q))-0.05, 1.05]);
title('Q fit: R^2 vs N');

%% ---- Residual histograms (pooled across all N) ----
figure('Color','w','Name','Residual histograms');
subplot(1,2,1);
histogram(detEtrue-detEmodel, 40);
grid on; xlabel('E_{true} - E_{model}'); ylabel('count');
title('Energy residuals, pooled over all N and cases');
subplot(1,2,2);
histogram(reshape(detQtrue-detQmodel,[],1), 60);
grid on; xlabel('Q_{true} - Q_{model}'); ylabel('count');
title('Q residuals (all 6 components), pooled over all N and cases');

end