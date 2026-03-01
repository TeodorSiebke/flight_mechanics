% calculate_total_np_nll.m
% Calculates the total aircraft Neutral Point using NLL-derived surface properties
% and Flow5-derived efficiency/downwash constants.

clear all;
clc;

fprintf('=== WINDEX TOTAL NEUTRAL POINT CALCULATION (NLL DATA) ===\n\n');

%% 1. INPUT DATA (from NLL Analysis)
% Wing properties (from analyze_wing.m)
S_w       = 7.2537;      % m^2
mac_w     = 0.6274;      % m
xle_mac_w = 0.0468;      % m (from Root LE)
cla_w     = 6.0411;      % /rad
h_nwb     = 0.2737;      % AC position as fraction of MAC_w

% Tail properties (from analyze_tail.m)
S_t       = 0.8600;      % m^2
mac_t     = 0.4149;      % m
xle_mac_t = 2.8166;      % m (from Root LE wing)
cla_t     = 4.9535;      % /rad
h_nt      = 0.2664;      % AC position as fraction of MAC_t

%% 2. AERO CONSTANTS (from windex_config.m)
eta       = 0.95;        % Tail efficiency
deps_da   = 0.176;       % Downwash gradient

%% 3. DERIVED GEOMETRY
% Absolute AC positions from Wing Root LE (Datum x=0)
Xac_w = xle_mac_w + h_nwb * mac_w;
Xac_t = xle_mac_t + h_nt * mac_t;

% Lever arm (Wing AC to Tail AC)
l_t = Xac_t - Xac_w;

% Tail Volume Coefficient
V_bar = (S_t * l_t) / (S_w * mac_w);

%% 4. STABILITY CALCULATIONS
% Total Aircraft Lift Slope
% CL_a = CL_a_w + (S_t/S_w) * eta * CL_a_t * (1 - deps_da)


cla_total = cla_w + (S_t/S_w) * eta * cla_t * (1-deps_da);

% Neutral Point (h_n) as fraction of MAC_w (referenced to Wing MAC LE)
% h_n = h_nwb + (1/CL_a_total) * V_bar * CL_a_t * eta * (1 - deps_da)
h_n_fraction = h_nwb + (1/cla_total) * V_bar * cla_t * eta * (1-deps_da);

% Neutral Point in meters from Wing Root LE
X_np = xle_mac_w + h_n_fraction * mac_w;

%% 5. OUTPUT RESULTS
fprintf('Wing AC (X_ac_w):    %.4f m (from Root LE)\n', Xac_w);
fprintf('Tail AC (X_ac_t):    %.4f m (from Root LE)\n', Xac_t);
fprintf('Tail Lever Arm (lt): %.4f m\n', l_t);
fprintf('Tail Volume (Vbar):  %.4f\n\n', V_bar);

fprintf('Combined CL_alpha:   %.4f /rad\n', cla_total);
fprintf('----------------------------------------------------\n');
fprintf('Neutral Point (h_n): %.4f (%% of MAC: %.2f%%)\n', h_n_fraction, h_n_fraction*100);
fprintf('Neutral Point (X_np): %.4f m (from Wing Root LE)\n', X_np);
fprintf('----------------------------------------------------\n');

% Comparison with h_nwb
fprintf('Static Margin (if CG at AC_w): %.2f%% MAC\n', (h_n_fraction - h_nwb)*100);
