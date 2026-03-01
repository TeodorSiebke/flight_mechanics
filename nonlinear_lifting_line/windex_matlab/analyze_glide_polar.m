% analyze_glide_polar.m
% Calculates and plots the trimmed glide polar (Sink Rate vs Airspeed)
% for different flap settings: -10, 0, 10, 20 degrees.

clear all; close all; clc;

% Add path to NLL library
addpath('../matlab');

% --- Setup ---
fprintf('Initializing Windex Model...\n');
acg = make_windex();

% Constants
mass = 300; % kg
g = 9.81;
W = mass * g;
rho = 1.225; % kg/m^3 (SL)
S_ref = acg.Sref;
c_ref = acg.cref;

% Find Neutral Point once (reference for CG)
fprintf('Calculating Neutral Point...\n');
[h_n, X_LE_wing_mac] = calculate_nll_np_total(acg);

% Set CG for a fixed Static Margin (e.g., 9%)
sm_target = 0.09;
h_cg = h_n - sm_target;
X_cg = X_LE_wing_mac + h_cg * c_ref;
acg.pref = [X_cg, 0, 0];
fprintf('Using SM = %.1f%% (CG at %.4f m)\n', sm_target*100, X_cg);

% --- Analysis Matrix ---
flap_vec = [0, 10, 20]; % [deg]
V_vec = 18:2:65;             % [m/s]
n_flap = length(flap_vec);
n_v = length(V_vec);

results = struct();

% Optimization options
options = optimset('Display','off', 'TolX', 1e-4, 'MaxIter', 100);

fprintf('\nStarting Trimmed Glide Polar Analysis...\n');

for k = 1:n_flap
    flap_deg = flap_vec(k);
    fprintf('  Flaps: %3d deg... ', flap_deg);
    
    res_V = [];
    res_Vs = [];
    res_al = [];
    res_de = [];
    res_LD = [];
    
    % Initial guesses
    al_guess = 2 * pi/180; 
    de_guess = 0;
    
    for i = 1:n_v
        V = V_vec(i);
        q = 0.5 * rho * V^2;
        CL_req = W / (q * S_ref);
        
        % Check if CL is reasonable before solving
        if CL_req > 2.0 || CL_req < -0.5
            continue; 
        end
        
        % Solve for [alpha, delta_e]
        % Unknowns x = [alpha_rad, delta_e_rad]
        obj_fun = @(x) trim_solver(x, acg, V, CL_req, flap_deg);
        
        try
            [x_sol, fval, exitflag] = fsolve(obj_fun, [al_guess, de_guess], options);
        catch
            exitflag = -1;
        end
        
        if exitflag > 0
            al_rad = x_sol(1);
            de_rad = x_sol(2);
            
            % Get final coefficients to calculate CD
            fs = set_fs_glide(acg, V, al_rad, de_rad, flap_deg);
            [CL, CD, ~] = plainpolar(acg, fs, al_rad);
            
            % Sink speed Vs = V * sin(gamma) approx V * (D/L) = V * (CD/CL)
            Vs = V * (CD(1) / CL(1));
            
            res_V(end+1) = V;
            res_Vs(end+1) = Vs;
            res_al(end+1) = al_rad * 180/pi;
            res_de(end+1) = de_rad * 180/pi;
            res_LD(end+1) = CL(1) / CD(1);
            
            % Update guesses
            al_guess = al_rad;
            de_guess = de_rad;
        end
    end
    
    results(k).flap = flap_deg;
    results(k).V = res_V;
    results(k).Vs = res_Vs;
    results(k).LD = res_LD;
    fprintf('Done (%d speeds).\n', length(res_V));
end

% --- Plotting ---
figure('Name', 'Trimmed Glide Polar', 'Color', 'w', 'Position', [100 100 800 600]);
colors = {'b', 'k', 'r', 'g'};
markers = {'o', 's', '^', 'd'};

hold on; grid on; box on;
all_V = [];
all_Vs = [];

for k = 1:n_flap
    if isempty(results(k).V), continue; end
    
    % Note: Sink rate is usually plotted positive downward or negative upward.
    % We will plot Vs on Y axis (positive = down) and V on X axis.
    plot(results(k).V, results(k).Vs, [colors{k} markers{k} '-'], 'LineWidth', 1.5, ...
        'DisplayName', ['Flaps ' num2str(results(k).flap) '^\circ']);
    
    all_V = [all_V, results(k).V];
    all_Vs = [all_Vs, results(k).Vs];
end

% --- Find GLOBAL Best Glide (Across all Flaps) ---
max_LD_global = 0;
best_flap_idx = 0;
best_v_idx = 0;

for k = 1:n_flap
    if isempty(results(k).LD), continue; end
    [local_max, local_idx] = max(results(k).LD);
    if local_max > max_LD_global
        max_LD_global = local_max;
        best_flap_idx = k;
        best_v_idx = local_idx;
    end
end

if best_flap_idx > 0
    V_best = results(best_flap_idx).V(best_v_idx);
    Vs_best = results(best_flap_idx).Vs(best_v_idx);
    flap_best = results(best_flap_idx).flap;
    
    % Draw tangent from origin (0,0) through (V_best, Vs_best)
    V_max_plot = max(V_vec)*1.1;
    V_ext = [0, V_max_plot];
    Vs_ext = V_ext * (Vs_best / V_best);
    plot(V_ext, Vs_ext, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    
    % Mark the point
    plot(V_best, Vs_best, 'ko', 'MarkerSize', 12, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Best Glide (L/D=%.1f at %d^\\circ flaps)', max_LD_global, flap_best));
    
    fprintf('\nGLOBAL Best Glide:\n');
    fprintf('  Flaps: %d deg\n', flap_best);
    fprintf('  Speed: %.1f m/s\n', V_best);
    fprintf('  Sink Rate: %.2f m/s\n', Vs_best);
    fprintf('  L/D Max:   %.2f\n', max_LD_global);
end

set(gca, 'YDir', 'reverse'); % Traditional glide polar has sink rate increasing downwards
xlabel('Airspeed V [m/s]');
ylabel('Sink Rate V_s [m/s] (Positive Down)');
title(['Trimmed Glide Polar (Mass = ' num2str(mass) ' kg)']);
legend('Location', 'best');

if ~isempty(all_V)
    xlim([0, max(all_V)*1.1]);
    ylim([0, max(all_Vs)*1.2]);
end

saveas(gcf, 'Windex_Trimmed_Glide_Polar.png');

% Plot L/D as well
figure('Name', 'Glide Ratio', 'Color', 'w', 'Position', [950 100 800 600]);
hold on; grid on;
for k = 1:n_flap
    if isempty(results(k).V), continue; end
    plot(results(k).V, results(k).LD, [colors{k} markers{k} '-'], 'LineWidth', 1.5, ...
        'DisplayName', ['Flaps ' num2str(results(k).flap) '^\circ']);
end
xlabel('Airspeed V [m/s]');
ylabel('Glide Ratio L/D');
title(['Glide Ratio (Trimmed, Mass = ' num2str(mass) ' kg)']);
legend('Location', 'best');
saveas(gcf, 'Windex_Trimmed_LD.png');

fprintf('\nAnalysis Finished. Plots saved as Windex_Trimmed_Glide_Polar.png and Windex_Trimmed_LD.png\n');

%% HELPER FUNCTIONS
function F = trim_solver(x, acg, V, CL_req, flap_deg)
    alpha = x(1);
    de = x(2);
    
    fs = set_fs_glide(acg, V, alpha, de, flap_deg);
    [vCL, ~, vCM] = plainpolar(acg, fs, alpha);
    
    F(1) = vCL(1) - CL_req;
    F(2) = vCM(1); % Cm = 0
end

function fs = set_fs_glide(acg, V, alpha, de, flap_deg)
    fs = flight_state(acg, V, 0, 0, 0);
    fs.alpha = alpha;
    fs.delta_elevator = de;
    fs.delta_flap = flap_deg * pi/180;
    fs.delta_aileron = 0; fs.delta_rudder = 0;
end
