% compare_np.m
% Run this to compare NLLT and Analytical Neutral Point

addpath('nonlinear_lifting_line/windex_matlab');
addpath('nonlinear_lifting_line/windex_fuselagetest');
addpath('nonlinear_lifting_line/matlab');

fprintf('=== NEUTRAL POINT COMPARISON ===\n');

% 1. NLLT RESULTS
fprintf('Running NLLT NP calculation...\n');
acg = make_windex();
[h_n_nllt, x_le_mac_nllt, x_np_nllt] = calculate_nll_np_total(acg);
% Note: calculate_nll_np_total uses acg.cref (0.641) for normalization.

% 2. ANALYZE COMPONENTS (Wing only vs Total)
% Set up for Wing-only to find h_ac_wb
acg_wing = acg;
acg_wing.vxhtp = []; % Remove tail
acg_wing.vxvtp = []; 
[h_n_wing, x_le_mac_wing, x_np_wing] = calculate_nll_np_total(acg_wing);

% 3. ANALYZE TAIL EFFECTIVENESS
% To find deps_da and a_h, we can run a tail-only analysis if needed, 
% or just look at the slope differences.

% 4. ANALYTICAL ESTIMATE (from windex_config)
% Manual re-calc based on windex_config.m
S_w = 7.41;
S_h = 0.87;
MAC = 0.6119;
l_t = 2.681;
aw = 0.1003 * 180/pi;
ah = 0.0774 * 180/pi;
eta = 0.95;
deps_da = 0.176;
h_nwb = 0.2486;

CL_a_total = aw + (S_h/S_w) * eta * ah * (1 - deps_da);
V_h = (S_h * l_t) / (S_w * MAC);
h_n_analytical = h_nwb + (1/CL_a_total) * V_h * ah * eta * (1 - deps_da);

fprintf('\nSummary Results:\n');
fprintf('Analytical h_n: %.4f (%% MAC, MAC=0.6119)\n', h_n_analytical);
fprintf('NLLT h_n:       %.4f (%% cref, cref=0.6410)\n', h_n_nllt);
fprintf('Wing-only h_n:  %.4f (This is h_ac_wb in NLLT)\n', h_n_wing);

% Check normalization factor
h_n_nllt_adj = (x_np_nllt - x_le_mac_nllt) / 0.6119;
fprintf('NLLT h_n (adj to MAC=0.6119): %.4f\n', h_n_nllt_adj);

% Calculate discrepancy
disc = h_n_nllt_adj - h_n_analytical;
fprintf('\nDiscrepancy (%% MAC): %.2f%%\n', disc * 100);

exit;
