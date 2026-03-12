% find_windex_np_datum.m
% Calculates the neutral point of the Windex from the datum (root LE).
% 1. Trim elevator for V=30 m/s
% 2. Generate CL and Cm polars
% 3. Linearize around alpha=0
% 4. Estimate NP from datum

clear all;
close all;
clc;

% Add paths for lifting-line functions (using absolute path for robustness)
addpath('c:\Users\teodo\Desktop\flight_mechanics\nonlinear_lifting_line\matlab');

fprintf('=== WINDEX NEUTRAL POINT ESTIMATION (FROM DATUM) ===\n\n');

%% 1. INITIALIZE GEOMETRY
fprintf('Initializing Windex geometry...\n');
acg = make_windex();

% Reference properties
S_ref = acg.Sref;
c_ref = acg.cref;
% Ensure reference point is at a realistic CG (e.g. 0.3m from Datum)
% to get a reasonable elevator deflection.
acg.pref = [0.34, 0, 0]; 
x_ref = acg.pref(1);

%% 2. TRIM ANALYSIS
uoo = 30; % m/s
mass = 300; % kg
g = 9.81;
rho = 1.225;
W = mass * g;
q = 0.5 * rho * uoo^2;
CL_req = W / (q * S_ref);

fprintf('Trimming for V = %.1f m/s (CL_req = %.4f)...\n', uoo, CL_req);

% Initial guesses [alpha_rad, delta_e_rad]
x0 = [2*pi/180, 0];
options = optimset('Display','off', 'TolX', 1e-6);

% Solve for trim [alpha, delta_e]
% We use the anonymous function defined below
obj_fun = @(x) trim_residuals(x, acg, uoo, CL_req);
try
    [x_sol, fval, exitflag] = fsolve(obj_fun, x0, options);
catch ME
    fprintf('Error during fsolve: %s\n', ME.message);
    rethrow(ME);
end

if exitflag > 0
    alpha_trim = x_sol(1);
    de_trim = x_sol(2);
    fprintf('Trim found:\n');
    fprintf('  Alpha:    %.3f deg\n', alpha_trim * 180/pi);
    fprintf('  Elevator: %.3f deg\n', de_trim * 180/pi);
else
    error('Trim failed to converge.');
end

%% 3. GENERATE POLARS
fprintf('\nGenerating polars around trim elevator setting...\n');

% Sweep alpha around 0 and trim
valpha_deg = -5:1:10;
valpha_rad = valpha_deg * pi/180;

% Create flight state with trimmed elevator
fs_trim = flight_state(acg, uoo, 0, 0, 0);
fs_trim.delta_elevator = de_trim;
fs_trim.delta_rudder = 0;
fs_trim.delta_flap = 0;
fs_trim.delta_aileron = 0;

[vCL, vCD, vCM] = plainpolar(acg, fs_trim, valpha_rad);

%% 4. LINEARIZE AROUND ALPHA = 2 deg
fprintf('Linearizing around alpha = 2 deg...\n');

% Linearize around alpha = 2 deg
alpha_lin_deg = 2.0;
[~, idx] = min(abs(valpha_deg - alpha_lin_deg));

% Use points +/- 1 deg around target
idx_m = idx - 1;
idx_p = idx + 1;

if idx_m < 1 || idx_p > length(valpha_rad)
    indices = [idx, idx+1];
else
    indices = [idx_m, idx_p];
end

CLa = (vCL(indices(2)) - vCL(indices(1))) / (valpha_rad(indices(2)) - valpha_rad(indices(1)));
CMa = (vCM(indices(2)) - vCM(indices(1))) / (valpha_rad(indices(2)) - valpha_rad(indices(1)));

fprintf('Slopes at alpha=%.1f deg:\n', alpha_lin_deg);
fprintf('  dCL/da: %.4f per rad (%.4f per deg)\n', CLa, CLa * pi/180);
fprintf('  dCm/da: %.4f per rad (%.4f per deg)\n', CMa, CMa * pi/180);

%% 5. ESTIMATE NEUTRAL POINT
% x_np = x_ref - c_ref * (dCm/da / dCL/da)
X_np = x_ref - c_ref * (CMa / CLa);

fprintf('\n=== RESULTS (Linearized at alpha=%.1f deg) ===\n', alpha_lin_deg);
fprintf('Neutral Point (X_np): %.4f m from Datum (Root LE)\n', X_np);
fprintf('Neutral Point (%% MAC): %.2f%%\n', (X_np / c_ref) * 100);

%% 6. PLOTTING
figure('Name', 'Windex Polars', 'Color', 'w', 'Position', [100 100 800 400]);

subplot(1,2,1);
plot(valpha_deg, vCL, 'b-o', 'LineWidth', 1.5);
hold on;
yline(CL_req, 'r--', 'Trim CL');
xline(alpha_lin_deg, 'k--', 'Lin Alpha');
xlabel('\alpha [deg]'); ylabel('C_L');
title(['Lift Polar (Lin at \alpha=' num2str(alpha_lin_deg) ')']);
grid on;

subplot(1,2,2);
plot(valpha_deg, vCM, 'm-s', 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'Cm=0');
xline(alpha_trim*180/pi, 'g--', 'Trim \alpha');
xline(alpha_lin_deg, 'k--', 'Lin Alpha');
xlabel('\alpha [deg]'); ylabel('C_m');
title(['Mom. Polar about CG=' num2str(x_ref) ' (Lin at \alpha=' num2str(alpha_lin_deg) ')']);
grid on;

saveas(gcf, 'Windex_NP_Polars.png');
fprintf('\nPolars saved as Windex_NP_Polars.png\n');


%% HELPER FUNCTIONS
function F = trim_residuals(x, acg, V, CL_req)
    alpha = x(1);
    de = x(2);
    
    fs = flight_state(acg, V, 0, alpha, 0);
    fs.delta_elevator = de;
    fs.delta_rudder = 0; fs.delta_flap = 0; fs.delta_aileron = 0;
    
    % Get coefficients at this state
    % plainpolar expects a flight condition and a vector of alphas
    % but we just want one point
    [CL, ~, CM] = plainpolar(acg, fs, alpha);
    
    F(1) = CL - CL_req;
    F(2) = CM; % Pitch trim about the current acg.pref (Datum)
end
