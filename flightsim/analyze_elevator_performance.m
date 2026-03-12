% analyze_elevator_performance.m
% Linearization and Performance Analysis of the Elevator PID Controller
% 
% This script investigates the stability and performance of the 
% altitude-hold PID controller across a range of altitudes and speeds.

clear; close all; clc;
fsm = make_fsim();

% --- 1. SETTINGS & ENVELOPE ---
v_vec = 25:5:100;       % Airspeed range (m/s)
h_vec = 500:500:10000;  % Altitude range (m)
n_v = length(v_vec);
n_h = length(h_vec);

% Controller Gains (Matches elevator_level_flight.m)
% Optimized for 30 m/s stability at 1000m (T2D > 200s)
Kp_h = 0.005093;
Ki_h = 0.000209;
Kd_h = 0.011718;
Kq   = 0.418435;

% Storage for metrics
% 1. Damping ratio of the dominant mode
% 2. Time-to-half (s) for stability visualization
damping_matrix = zeros(n_h, n_v);
t_half_matrix = zeros(n_h, n_v);
stability_flag = zeros(n_h, n_v);

fprintf('Starting Envelope Sweep (%d points)...\n', n_v * n_h);

% --- 2. ENVELOPE SWEEP ---
for i = 1:n_h
    h_target = h_vec(i);
    for j = 1:n_v
        v_target = v_vec(j);
        
        % A. TRIM THE AIRCRAFT
        % [u w q theta x h fuel de dp]
        x_tr = [v_target 0.5 0 0.05 0 h_target 0 0.06 0.1]';
        ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8];
        xt = x_tr(ivar);
        converged = false;
        for iter = 1:50
            x_tr(ivar) = xt;
            xd = fplmod(0, x_tr, fsm);
            res = xd(ifun); res(5) = res(5) - v_target;
            if norm(res) < 1e-10, converged = true; break; end
            xs = 1e-7; J = zeros(5,5);
            for k = 1:5
                xh = x_tr; xh(ivar(k)) = xh(ivar(k)) + xs;
                xmh = x_tr; xmh(ivar(k)) = xmh(ivar(k)) - xs;
                xdh = fplmod(0, xh, fsm); xdm = fplmod(0, xmh, fsm);
                J(:,k) = (xdh(ifun) - xdm(ifun)) / (2*xs);
            end
            xt = xt - 0.7 * (J \ res);
        end
        x_tr(ivar) = xt;
        
        if ~converged
            fprintf('  V=%.1f H=%.0f: FAILED TO TRIM\n', v_target, h_target);
            stability_flag(i,j) = 0; % Mark as unstable
            damping_matrix(i,j) = 0; % No damping if not trimmed
            continue;
        end
        
        % B. LINEARIZATION
        % Augmented State: [u w q theta h h_int]
        % (We ignore x_dist and fuel for stability analysis)
        % x_aug = [x(1) x(2) x(3) x(4) x(6) x_int]
        state_idx = [1 2 3 4 6];
        n_aug = 6;
        A = zeros(n_aug, n_aug);
        xs_lin = 1e-5;
        
        % Construct the wrapper for the closed-loop system
        % x_vec: [u w q theta h h_int]
        cl_system = @(x_vec) closed_loop_dynamics(x_vec, x_tr, fsm, Kp_h, Ki_h, Kd_h, Kq, h_target);
        
        for k = 1:n_aug
            xh = [x_tr(state_idx); 0]; xh(k) = xh(k) + xs_lin;
            xmh = [x_tr(state_idx); 0]; xmh(k) = xmh(k) - xs_lin;
            
            xdh = cl_system(xh);
            xdm = cl_system(xmh);
            A(:,k) = (xdh - xdm) / (2*xs_lin);
        end
        
        % C. PERFORMANCE ANALYSIS
        eg = eig(A);
        real_parts = real(eg);
        imag_parts = imag(eg);
        
        % Soft Stability: Accept real parts if Time-to-Double > 200s
        % Threshold = ln(2)/200 approx 0.00346
        stability_threshold = log(2)/200;
        stability_flag(i,j) = all(real_parts < stability_threshold); 
        
        % Metric 1: Minimum Damping Ratio (Oscillatory modes)
        oscillatory = find(abs(imag_parts) > 1e-4); % Catch slow Phugoid
        if ~isempty(oscillatory)
            zetas = -real_parts(oscillatory) ./ abs(eg(oscillatory));
            damping_matrix(i,j) = min(zetas);
        else
            damping_matrix(i,j) = 1.0; % Overdamped
        end
        
        % Metric 2: Growth/Decay Rate (1/Time value)
        % We use lambda / ln(2) which is 1/T
        [~, dom_idx] = max(real_parts);
        t_half_matrix(i,j) = real_parts(dom_idx) / log(2); 
        
        fprintf('  V=%.1f H=%.0f: Stable=%d, Min Damping=%.3f, 1/T=%.4f\n', ...
            v_target, h_target, stability_flag(i,j), damping_matrix(i,j), t_half_matrix(i,j));
    end
end

% --- 3. VISUALIZATION ---
figure('Name', 'Control System Performance Sweep', 'Color', 'w', 'Position', [100 100 1100 500]);
[V_grid, H_grid] = meshgrid(v_vec, h_vec);

subplot(1,2,1);
imagesc(v_vec, h_vec, damping_matrix);
set(gca, 'YDir', 'normal');
colorbar; colormap(gca, parula);
title('Damping Ratio \zeta (Higher is Better)');
xlabel('Airspeed (m/s)'); ylabel('Altitude (m)');
clim([0 1]); % Standard damping range
xticks(v_vec); % Show all airspeeds on X-axis

% Add Acceptable Envelope Contour (e.g., Zeta >= 0.3)
hold on;
[C, h_cont] = contour(V_grid, H_grid, damping_matrix, [0.3 0.3], 'r', 'LineWidth', 2.5);
clabel(C, h_cont, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
legend(h_cont, 'Acceptable Envelope (\zeta \geq 0.3)', 'Location', 'southoutside');

subplot(1,2,2);
% Plot Stability Intensity Map (1/T)
% Stability logic:
% Stable (lambda < 0): Display as -1/T_half -> highly negative is Deep Green
% Unstable (lambda > 0): Display as 1/T_double -> highly positive is Deep Red
% Neutral (lambda = 0): Display as 0 -> White
intensity_map = t_half_matrix;
imagesc(v_vec, h_vec, intensity_map);
set(gca, 'YDir', 'normal');

% Custom Divergent Colormap: Green (Stable) -> White (Neutral) -> Red (Unstable)
% We want highly negative (Green) -> zero (White) -> highly positive (Red)
n_colors = 128;
cmap_g = [linspace(0,1,n_colors/2)', linspace(1,1,n_colors/2)', linspace(0,1,n_colors/2)']; % Green to White
cmap_r = [linspace(1,1,n_colors/2)', linspace(1,0,n_colors/2)', linspace(1,0,n_colors/2)']; % White to Red
cmap = [cmap_g; cmap_r];
colormap(gca, cmap);

% Symmetric limit to center white at 0 (Inf seconds)
% Set limit to 0.1 (Rate = 1/10s) so fast modes (>10s) are saturated
% and slow modes (100s = 0.01) are faint (10% intensity)
max_viz = 0.1; 
clim([-max_viz max_viz]);

cb = colorbar;
% Custom labels for the colorbar to show seconds
ticks = [-0.1 -0.05 -0.01 0 0.01 0.05 0.1];
tick_labels = cell(size(ticks));
for k = 1:length(ticks)
    if ticks(k) < -0.0001
        tick_labels{k} = sprintf('T_{1/2}=%.0fs', 1/abs(ticks(k)));
    elseif ticks(k) > 0.0001
        tick_labels{k} = sprintf('T_{2}=%.0fs', 1/ticks(k));
    else
        tick_labels{k} = 'Inf (Neutral)';
    end
end
set(cb, 'Ticks', ticks, 'TickLabels', tick_labels);

title('Stability (Fast = Deep Color, 100s+ = Faint)');
xlabel('Airspeed (m/s)'); ylabel('Altitude (m)');
xticks(v_vec); % Show all airspeeds on X-axis

sgtitle('Elevator PID Performance across Altitude-Speed Envelope');
saveas(gcf, 'control_performance_envelope.png');

fprintf('\nAnalysis complete. Heatmaps saved.\n');

% --- HELPER: CLOSED LOOP DYNAMICS ---
function xdot_aug = closed_loop_dynamics(x_vec, x_tr, fsm, Kp, Ki, Kd, Kq, h_target)
    % x_vec: [u w q theta h h_int]
    u = x_vec(1);
    w = x_vec(2);
    q = x_vec(3);
    theta = x_vec(4);
    h = x_vec(5);
    h_int = x_vec(6);
    
    de_trim = x_tr(8);
    dp_trim = x_tr(9);
    
    % PID Logic
    h_err = h_target - h;
    h_dot = - ( -sin(theta)*u + cos(theta)*w );
    de = de_trim - (Kp * h_err + Ki * h_int - Kd * h_dot) + Kq * q;
    
    % Full state for fplmod
    xf = [u; w; q; theta; 0; h; 0; de; dp_trim];
    xd_phys = fplmod(0, xf, fsm);
    
    % Augmented derivative
    xdot_aug = [xd_phys(1); xd_phys(2); xd_phys(3); xd_phys(4); xd_phys(6); h_err];
end
