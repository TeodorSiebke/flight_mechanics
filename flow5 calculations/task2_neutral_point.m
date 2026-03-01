%% Task 2: Neutral Point Calculation
% Calculates the analytical Neutral Point of the aircraft.
clear; clc;

% Load Configuration
if ~exist('Config', 'var')
    windex_config;
end

fprintf('\n=== TASK 2: NEUTRAL POINT CALCULATION ===\n');

% UNPACK CONFIG
% Making variables easier to read for the equations
A = Config.Aero;

% Derived values
V_bar_ht = (A.S_tail * A.l_t) / (A.S_wing * A.MAC); % Horizontal tail volume coefficient

% CALCULATIONS
% Total Aircraft Lift-Curve Slope (Eq 3)
% CL_alpha = CL_alpha_wb + (S_ht/S) * eta_ht * CL_alpha_ht * (1 - deps_dalpha)
CL_alpha = A.CL_alpha_wb + (A.S_tail/A.S_wing) * A.eta_ht * A.CL_alpha_ht * (1 - A.deps_dalpha);

% Neutral Point Calculation (h_n)
% h_n = h_nwb + (1/CL_alpha) * V_bar_ht * CL_alpha_ht * eta_ht * (1 - deps_dalpha)
h_n = A.h_nwb + (1/CL_alpha) * V_bar_ht * A.CL_alpha_ht * A.eta_ht * (1 - A.deps_dalpha);

% DISPLAY RESULTS
fprintf('----------------------------------------------------\n');
fprintf('Wing Lift Slope (CL_a_wb):   %.4f /rad\n', A.CL_alpha_wb);
fprintf('Tail Volume Coeff (V_bar):   %.4f\n', V_bar_ht);
fprintf('----------------------------------------------------\n');
fprintf('Total Lift-Curve Slope:      %.4f /rad\n', CL_alpha);
fprintf('Neutral Point (h_n):         %.4f (fraction of MAC)\n', h_n);
fprintf('----------------------------------------------------\n');
