% analyze_trim_nll.m
% Trimmed Lift Curve Analysis using Nonlinear Lifting Line Theory (NLLT).
% Revamped to use full geometry from make_windex.m and perform nonlinear trim.

clear all;
close all;
clc;

addpath('../matlab');

fprintf('=== WINDEX TRIM ANALYSIS (FULL NONLINEAR) ===\n\n');

%% 1. MODEL SETUP & GEOMETRY
fprintf('Initializing Windex geometry...\n');
acg = make_windex();

% Reference properties from make_windex (User request)
S_ref = acg.Sref;
c_ref = acg.cref; 
b_ref = acg.bref;

% Numerical Wing MAC LE (needed for X_ac reference)
rw = acg.vxwing;
l_seg = sqrt(sum((acg.pb - acg.pa).^2, 2));
S_segs = acg.chord .* l_seg;
S_wing = sum(S_segs(rw));
X_LE_wing_mac = sum(acg.pa(rw,1) .* S_segs(rw)) / S_wing;

fprintf('Reference Geometry:\n');
fprintf('S_ref: %.4f m^2\n', S_ref);
fprintf('c_ref: %.4f m\n', c_ref);
fprintf('MAC LE (calc): %.4f m\n', X_LE_wing_mac);

%% 2. CALCULATE NEUTRAL POINT (NP)
% Find NP by finding the reference point that makes dCm/dCL = 0 in the linear range.
fprintf('\nCalculating Neutral Point (NP)...\n');

% Set temporary reference to Wing MAC LE
acg.pref = [X_LE_wing_mac, 0, 0]; 

uoo = 30;
alpha_np = [-1, 1] * pi/180; % Small linear range
fs0 = flight_state(acg, uoo, 0, 0, 0);
% Zeros for all controls
fs0.delta_elevator = 0; fs0.delta_rudder = 0; fs0.delta_flap = 0; fs0.delta_aileron = 0;

[CL_np, ~, CM_np] = plainpolar(acg, fs0, alpha_np);

cla_np = (CL_np(2) - CL_np(1)) / (alpha_np(2) - alpha_np(1));
cma_np = (CM_np(2) - CM_np(1)) / (alpha_np(2) - alpha_np(1));
dcm_dcl = cma_np / cla_np;

% X_np = X_ref - (dCm/dCL) * c_ref
X_np = X_LE_wing_mac - dcm_dcl * c_ref;
h_n = (X_np - X_LE_wing_mac) / c_ref;

fprintf('Neutral Point (X_np): %.4f m (from Root LE)\n', X_np);
fprintf('Neutral Point (h_n):  %.4f (%% MAC)\n', h_n);

%% 3. SET CG TARGET
% User requested 9% Static Margin? Let's use 10% as a robust default or stick to 9% if previously used.
sm_target = 0.09; 
h_cg = h_n - sm_target;
X_cg = X_LE_wing_mac + h_cg * c_ref;

fprintf('Target Static Margin: %.1f%%\n', sm_target * 100);
fprintf('Target CG (X_cg):     %.4f m\n', X_cg);

% UPDATE REFERENCE POINT TO CG
acg.pref = [X_cg, 0, 0];

%% 4. TRIM VS AIRSPEED ANALYSIS
% Loop over Static Margins and Airspeeds
mass = 300; % kg
g = 9.81;
W = mass * g;
rho = 1.225; % Sea level density

sm_vec = [0.03, 0.09, 0.15, 0.30];
V_vec = [22:0.5:24, 26:2:60]; % m/s
n_sm = length(sm_vec);
n_v = length(V_vec);

% Storage
trim_results = struct();

fprintf('\nStarting Trim vs Airspeed Analysis (Mass = %.0f kg)...\n', mass);

for k = 1:n_sm
    current_sm = sm_vec(k);
    h_cg = h_n - current_sm;
    X_cg = X_LE_wing_mac + h_cg * c_ref;
    
    % Update CG
    acg.pref = [X_cg, 0, 0];
    
    fprintf('\nConfig %d: SM = %.0f%%, CG = %.4f m\n', k, current_sm*100, X_cg);
    fprintf('%-10s | %-10s | %-10s | %-10s\n', 'V [m/s]', 'CL_req', 'Alpha', 'Delta_e');
    
    res_V = zeros(1, n_v);
    res_de = zeros(1, n_v);
    res_al = zeros(1, n_v);
    
    % Initial guesses
    al_guess = 0; de_guess = 0;
    % Limit iterations to prevent hanging
    options = optimset('Display','off', 'TolX', 1e-4, 'MaxIter', 50, 'MaxFunEvals', 200);
    
    for i = 1:n_v
        V = V_vec(i);
        q = 0.5 * rho * V^2;
        CL_req = W / (q * S_ref);
        
        fprintf('  Solving for V = %4.1f m/s (CL_req = %.4f)... ', V, CL_req);
        
        % Solve for [alpha, delta_e]
        % Unknowns x = [alpha_rad, delta_e_rad]
        
        obj_fun = @(x) trim_solver(x, acg, V, CL_req);
        
        try
            [x_sol, fval, exitflag] = fsolve(obj_fun, [al_guess, de_guess*pi/180], options);
        catch
            exitflag = -1;
        end
        
        if exitflag > 0
            al_rad = x_sol(1);
            de_rad = x_sol(2);
            
            res_V(i) = V;
            res_al(i) = al_rad * 180/pi;
            res_de(i) = de_rad * 180/pi;
            
            % Update guess
            al_guess = al_rad;
            de_guess = res_de(i);
            
            fprintf('Done. Alpha=%.2f, De=%.2f\n', res_al(i), res_de(i));
        else
             fprintf('Failed or did not converge.\n');
             res_V(i) = V; res_al(i)=NaN; res_de(i)=NaN;
        end
    end
    
    trim_results(k).SM = current_sm;
    trim_results(k).V = res_V;
    trim_results(k).de = res_de;
    trim_results(k).alpha = res_al;
end

%% 5. RE-RUN ORIGINAL TRIM CURVE FOR PLOT 1 (Using SM=9% as reference)
% We need the Trimmed CL vs Alpha curve for the first plot requested previously.
% We can just use the results from the SM=9% case (index 2) or re-run a detailed alpha sweep.
% The user wants "Trimmed Lift Slope vs Alpha".
% Let's re-run the detailed alpha sweep for SM=0.09 for consistency with previous request.

sm_ref = 0.09;
h_cg_ref = h_n - sm_ref;
X_cg_ref = X_LE_wing_mac + h_cg_ref * c_ref;
acg.pref = [X_cg_ref, 0, 0];

alpha_deg_sweep = -4:1:12;
alpha_rad_sweep = alpha_deg_sweep * pi/180;
n_pts = length(alpha_deg_sweep);
trim_CL_ref = zeros(1, n_pts);
uoo = 30; % Reference speed for CL curve

% Untrimmed
fs_un = flight_state(acg, uoo, 0, 0, 0);
fs_un.delta_elevator = 0; fs_un.delta_rudder = 0; fs_un.delta_flap = 0; fs_un.delta_aileron = 0;
[CL_un, ~, ~] = plainpolar(acg, fs_un, alpha_rad_sweep);

% Trimmed
de_guess = 0;
for i = 1:n_pts
    al = alpha_rad_sweep(i);
    obj_fun_cm = @(de_rad) get_cm_val(acg, uoo, al, de_rad);
    [de_sol_rad, ~, exitflag] = fsolve(obj_fun_cm, de_guess*pi/180, options);
    [vCL, ~, ~] = plainpolar(acg, set_fs(acg, uoo, al, de_sol_rad), al);
    trim_CL_ref(i) = vCL(1);
    de_guess = de_sol_rad * 180/pi;
end
dCL_da_trim = gradient(trim_CL_ref, alpha_rad_sweep);


%% 6. PLOTTING

% Figure 1: Trimmed Lift Curve (Ref SM=9%)
figure('Name', 'Trimmed Lift Curve', 'Color', 'w', 'Position', [100 100 600 500]);
plot(alpha_deg_sweep, trim_CL_ref, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Trimmed (Cm=0)');
hold on;
plot(alpha_deg_sweep, CL_un, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Untrimmed (\delta_e=0)');
grid on; xlabel('\alpha [deg]'); ylabel('C_L');
title(['Lift Curve (U = ' num2str(uoo) ' m/s, SM = ' num2str(sm_ref*100) '%)']);
legend('Location', 'best');
saveas(gcf, 'Windex_Trim_CL_NLL.png');

% Figure 2: Elevator Deflection vs Airspeed (Multi-SM)
figure('Name', 'Trim Elevator vs Speed', 'Color', 'w', 'Position', [750 100 600 500]);
colors = ['r', 'b', 'g', 'm'];
hold on;
for k = 1:n_sm
    plot(trim_results(k).V, trim_results(k).de, [colors(k) '-o'], 'LineWidth', 1.5, ...
        'DisplayName', ['SM = ' num2str(trim_results(k).SM*100) '%']);
end
grid on; xlabel('Airspeed [m/s]'); ylabel('Elevator Deflection [deg]');
title(['Trim Elevator for Mass = ' num2str(mass) ' kg']);
yline(10, 'k--', 'DisplayName', '\delta_e = 10 data'); 
yline(-10, 'k--', 'DisplayName', '\delta_e = -10 data');
legend('Location', 'best');
saveas(gcf, 'Windex_Trim_DeltaE_vs_Speed.png');

fprintf('\nDone. Figures saved.\n');

%% HELPER FUNCTIONS
function F = trim_solver(x, acg, V, CL_req)
    alpha = x(1);
    de = x(2);
    
    fs = set_fs(acg, V, alpha, de);
    [vCL, vCD, vCM] = plainpolar(acg, fs, alpha);
    
    F(1) = vCL(1) - CL_req;
    F(2) = vCM(1); % Cm = 0
end

function cm = get_cm_val(acg, V, alpha, de)
    fs = set_fs(acg, V, alpha, de);
    [~, ~, vCM] = plainpolar(acg, fs, alpha);
    cm = vCM(1);
end

function fs = set_fs(acg, V, alpha, de)
    fs = flight_state(acg, V, 0, 0, 0);
    fs.alpha = alpha;
    fs.delta_elevator = de;
    fs.delta_rudder = 0; fs.delta_flap = 0; fs.delta_aileron = 0;
end
