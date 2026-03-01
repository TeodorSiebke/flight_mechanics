% compare_np_detailed.m
% Detailed comparison of NLLT components vs Analytical assumptions

addpath('../matlab');
addpath('../../flow5 calculations');

fprintf('\n=== DETAILED NEUTRAL POINT DIAGNOSTICS ===\n');

% 1. GEOMETRY CHECK
acg = make_windex();
windex_config; % Loads Config
A = Config.Aero;

fprintf('\nGEOMETRY:\n');
fprintf('  NLLT S_ref: %.4f, cref: %.4f\n', acg.Sref, acg.cref);
fprintf('  Analytical S: %.4f, MAC: %.4f\n', A.S_wing, A.MAC);

% 2. ALPHA SWEEP - WING ONLY
fprintf('\nRUNNING NLLT (WING ONLY)...\n');
acg_wing = acg;
acg_wing.vxhtp = []; % Remove tail
acg_wing.vxvtp = []; 

uoo = 30;
alpha_rad = [-1, 1] * pi/180;
fs = flight_state(acg_wing, uoo, 0, 0, 0); 
[CL_w, ~, CM_w] = plainpolar(acg_wing, fs, alpha_rad);

cla_w = (CL_w(2) - CL_w(1)) / (alpha_rad(2) - alpha_rad(1));
cma_w = (CM_w(2) - CM_w(1)) / (alpha_rad(2) - alpha_rad(1));

% NLLT Wing AC (from Root LE, ref acg.pref = [0,0,0])
x_ac_w = - (cma_w / cla_w) * acg.cref; 
% User's MAC LE calculation logic from calculate_nll_np_total
rw = acg.vxwing;
l_seg = sqrt(sum((acg.pb - acg.pa).^2, 2));
S_segs = acg.chord .* l_seg;
S_wing = sum(S_segs(rw)); 
X_LE_wing_mac = sum(acg.pa(rw,1) .* S_segs(rw)) / S_wing;

h_ac_w = (x_ac_w - X_LE_wing_mac) / acg.cref;
h_ac_w_mac = (x_ac_w - X_LE_wing_mac) / A.MAC;

fprintf('  Wing Lift Slope (NLLT): %.4f /rad (vs Analytical %.4f)\n', cla_w, A.CL_alpha_wb);
fprintf('  Wing AC (NLLT): %.4f %% cref (%.4f %% MAC)\n', h_ac_w*100, h_ac_w_mac*100);
fprintf('  Wing AC (Analytical): %.1f %% MAC\n', A.h_nwb * 100);

% 3. ALPHA SWEEP - TOTAL
fprintf('\nRUNNING NLLT (TOTAL)...\n');
fs_tot = flight_state(acg, uoo, 0, 0, 0);
[CL_t, ~, CM_t] = plainpolar(acg, fs_tot, alpha_rad);

cla_t = (CL_t(2) - CL_t(1)) / (alpha_rad(2) - alpha_rad(1));
cma_t = (CM_t(2) - CM_t(1)) / (alpha_rad(2) - alpha_rad(1));

x_np_t = - (cma_t / cla_t) * acg.cref;
h_n_t = (x_np_t - X_LE_wing_mac) / acg.cref;
h_n_t_mac = (x_np_t - X_LE_wing_mac) / A.MAC;

fprintf('  Total Lift Slope (NLLT): %.4f /rad (vs Analytical %.4f)\n', cla_t, A.CL_alpha_wb + (A.S_tail/A.S_wing)*A.eta_ht*A.CL_alpha_ht*(1-A.deps_dalpha));
fprintf('  Total NP (NLLT): %.4f %% cref (%.4f %% MAC)\n', h_n_t*100, h_n_t_mac*100);

% 4. TAIL CONTRIBUTION ANALYSIS (Estimate)
% dCL_tail = CL_total - CL_wing (at alpha=0, approx)
% But let's look at slope difference
cla_h_eff = (cla_t - cla_w) / (A.S_tail / A.S_wing); % This is eta * a_h * (1 - deps_da)
fprintf('\nTAIL EFFECTIVENESS:\n');
fprintf('  Effective Tail Slope [eta*ah*(1-deps)]: %.4f (NLLT) vs %.4f (Analytical)\n', ...
    cla_h_eff, A.eta_ht * A.CL_alpha_ht * (1 - A.deps_dalpha));

% Estimate Downwash
% Need Tail-only run to get ah
acg_tail = acg;
acg_tail.vxwing = []; 
acg_tail.vxvtp = [];
fs_tail = flight_state(acg_tail, uoo, 0, 0, 0);
[CL_h, ~, ~] = plainpolar(acg_tail, fs_tail, alpha_rad);
cla_h = (CL_h(2) - CL_h(1)) / (alpha_rad(2) - alpha_rad(1)) * (acg.Sref / A.S_tail);

fprintf('  Isolated Tail Slope (a_h): %.4f /rad (NLLT) vs %.4f (Analytical)\n', cla_h, A.CL_alpha_ht);

% One-minus-deps_da
one_minus_deps = cla_h_eff / (A.eta_ht * cla_h);
fprintf('  Estimated (1 - deps/da): %.4f (if eta=%.2f) vs Analytical %.4f\n', ...
    one_minus_deps, A.eta_ht, 1 - A.deps_dalpha);

fprintf('\nCONCLUSION:\n');
diff_h_n = h_n_t_mac - A.h_nwb - (1/cla_t) * (A.S_tail*A.l_t/(A.S_wing*A.MAC)) * cla_h * A.eta_ht * (1-A.deps_dalpha);
fprintf('  NLLT shifted NP back by %.2f%% MAC relative to analytical components.\n', (h_n_t_mac*100 - Config.Aero.h_nwb*100));

exit;
