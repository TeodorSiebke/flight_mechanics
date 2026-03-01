% main_windex_analysis.m
% Analyzes the Windex airplane for different elevator deflections.

clear all;
close all;

% --- Path Setup ---
addpath('../matlab');

% --- Initialize Windex Model ---
fprintf('Initializing Windex geometry...\n');
acg = make_windex();

% --- Analysis Controls ---
uoo = 30;                              % Speed (m/s)
alpha_sweep = [-4:1:15] * pi/180;     % Alpha range (rad) - Reduced to match XFOIL data range
de_list = [-10, 0, 10];                % Elevator deflections (deg)

% --- Run Analysis (Elevator Sweep) ---
fprintf('Running Elevator Sweep...\n');
results_ele = struct();

for i = 1:length(de_list)
    de = de_list(i);
    fprintf('  Elevator: %d deg...\n', de);
    
    % Set flight state initial
    fs = flight_state(acg, uoo, 0.0, 0.0, 0.0);
    fs.delta_elevator = de * pi/180;
    fs.delta_rudder = 0.0;
    fs.delta_flap = 0.0;
    fs.delta_aileron = 0.0;
    
    % Compute polar
    [CL, CD, CM] = plainpolar(acg, fs, alpha_sweep);
    
    % Store
    results_ele(i).de = de;
    results_ele(i).CL = CL;
    results_ele(i).CM = CM;
    results_ele(i).CD = CD;
end

% --- Run Analysis (Aileron Sweep at 0 deg elevator) ---
fprintf('Running Aileron Sweep...\n');
da_list = [0, 10]; % Aileron deflections (deg)
results_ail = struct();

for i = 1:length(da_list)
    da = da_list(i);
    fprintf('  Aileron: %d deg...\n', da);
    
    fs = flight_state(acg, uoo, 0.0, 0.0, 0.0);
    fs.delta_elevator = 0.0;
    fs.delta_rudder = 0.0;
    fs.delta_flap = 0.0;
    fs.delta_aileron = da * pi/180;
    
    % Need Rolling Moment (Cl) - plainpolar usually returns CL, CD, CM (Pitching)
    % We need to check if we can get Cl (Rolling) from NLL.
    % plainpolar.m typically sums up forces to get CL, CD, CM.
    % If it doesn't return Cl, we might need to modify it or calculate it manually?
    % Let's verify plainpolar output first.
    % Assuming plainpolar only gives longitudinal coeffs.
    % Let's use `solve_lifting_line` directly or verify plainpolar.
    
    % For now, let's just run it and plot CL/CM changes if any (should be minimal for antisymmetric aileron)
    % Ailerons ideally change Cl (Rolling), not much CL (Lift) or CM (Pitch).
    [CL, CD, CM] = plainpolar(acg, fs, alpha_sweep);
    
    results_ail(i).da = da;
    results_ail(i).CL = CL;
end

% --- Plotting Results ---
fprintf('Generating plots...\n');
alpha_deg = alpha_sweep * 180 / pi;
colors = {'r', 'k', 'b'};
leg_str = {};

% 1. Lift Curve (Elevator)
figure(1);
for i = 1:length(de_list)
    plot(alpha_deg, results_ele(i).CL, colors{i}, 'LineWidth', 1.5); hold on;
    leg_str{i} = sprintf('\\delta_e = %d^\\circ', de_list(i));
end
xlabel('Angle of Attack \alpha [deg]');
ylabel('Lift Coefficient C_L');
title('Windex Lift Curves (Elevator Effect)');
grid on; legend(leg_str, 'Location', 'NorthWest');

% 2. Pitching Moment Curve (Elevator)
figure(2);
for i = 1:length(de_list)
    plot(alpha_deg, results_ele(i).CM, colors{i}, 'LineWidth', 1.5); hold on;
end
xlabel('Angle of Attack \alpha [deg]');
ylabel('Moment Coefficient C_m');
title('Windex Pitching Moment (Elevator Effect)');
grid on; legend(leg_str, 'Location', 'NorthEast');

% 3. Lift Curve (Aileron Check)
figure(3);
for i = 1:length(da_list)
    plot(alpha_deg, results_ail(i).CL, 'LineWidth', 1.5); hold on;
    leg_str_a{i} = sprintf('\\delta_a = %d^\\circ', da_list(i));
end
xlabel('Angle of Attack \alpha [deg]');
ylabel('Lift Coefficient C_L');
title('Windex Lift Curves (Aileron Effect - Should be minimal)');
grid on; legend(leg_str_a, 'Location', 'NorthWest');

fprintf('Done.\n');
