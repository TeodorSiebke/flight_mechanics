% compare_np.m
% Run this to compare NLLT and Analytical Neutral Point
% Run from: c:\Users\teodo\Desktop\flight_mechanics\nonlinear_lifting_line\windex_matlab

addpath('../matlab');
addpath('../../flow5 calculations');

fprintf('=== NEUTRAL POINT COMPARISON ===\n');

% 1. NLLT RESULTS
fprintf('Running NLLT NP calculation...\n');
acg = make_windex();
[h_n_nllt, x_le_mac_nllt, x_np_nllt] = calculate_nll_np_total(acg);

% 2. ANALYZE COMPONENTS (Wing only vs Total)
acg_wing = acg;
acg_wing.vxhtp = []; % Remove tail
acg_wing.vxvtp = []; 
[h_n_wing, x_le_mac_wing, x_np_wing] = calculate_nll_np_total(acg_wing);

% 3. ANALYTICAL ESTIMATE (from windex_config)
windex_config; % Loads Config
A = Config.Aero;

S_w = A.S_wing;
S_h = A.S_tail;
MAC = A.MAC;
l_t = A.l_t;
aw = A.CL_alpha_wb;
ah = A.CL_alpha_ht;
eta = A.eta_ht;
deps_da = A.deps_dalpha;
h_nwb = A.h_nwb;

CL_a_total = aw + (S_h/S_w) * eta * ah * (1 - deps_da);
V_h = (S_h * l_t) / (S_w * MAC);
h_n_analytical = h_nwb + (1/CL_a_total) * V_h * ah * eta * (1 - deps_da);

fprintf('\nSummary Results:\n');
fprintf('Analytical h_n: %.4f (%% MAC, MAC=%.4f)\n', h_n_analytical, MAC);
fprintf('NLLT h_n:       %.4f (%% cref, cref=%.4f)\n', h_n_nllt, acg.cref);
fprintf('Wing-only h_n:  %.4f (This is h_ac_wb in NLLT)\n', h_n_wing);

% Check normalization factor
h_n_nllt_adj = (x_np_nllt - x_le_mac_nllt) / MAC;
fprintf('NLLT h_n (adj to MAC=%.4f): %.4f\n', MAC, h_n_nllt_adj);

% Calculate discrepancy
disc = h_n_nllt_adj - h_n_analytical;
fprintf('\nDiscrepancy (%% MAC): %.2f%%\n', disc * 100);

exit;
