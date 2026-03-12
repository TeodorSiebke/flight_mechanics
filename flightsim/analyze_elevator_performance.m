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

% Controller Gains (Cruise Smoothed)
Kp_v = 0.001785;
Ki_v = 0.000259;
Kd_v = 0.018147;
Kq   = 0.054794;

% Storage for metrics
% 1. Damping ratio of the dominant mode
% 2. Time-to-half (s) for stability visualization
damping_matrix = NaN(n_h, n_v);
t_half_matrix = NaN(n_h, n_v);
stability_flag = zeros(n_h, n_v);

% Initial guess for the first point
x_tr = [v_vec(1) 0.5 0 0.05 0 h_vec(1) 0 0.06 0.1]';

fprintf('Starting Envelope Sweep (%d points)...\n', n_v * n_h);

% --- 2. ENVELOPE SWEEP ---
for i = 1:n_h
    h_target = h_vec(i);
    for j = 1:n_v
        v_target = v_vec(j);
        
        % A. TRIM THE AIRCRAFT
        % Standard guess for EACH point to prevent error propagation
        x_tr = [v_target 0.5 0 0.05 0 h_target 0 0.06 0.1]';
        ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8];
        xt = x_tr(ivar);
        converged = false;
        for iter = 1:50
            x_tr(ivar) = xt;
            xd = fplmod(0, x_tr, fsm);
            res = [xd(1); xd(2); xd(3); xd(6); sqrt(x_tr(1)^2+x_tr(2)^2)-v_target];
            if norm(res) < 1e-10, converged = true; break; end
            xs = 1e-8; J = zeros(5,5);
            for k = 1:5
                xh = x_tr; xh(ivar(k)) = xh(ivar(k)) + xs;
                xmh = x_tr; xmh(ivar(k)) = xmh(ivar(k)) - xs;
                xdh = fplmod(0, xh, fsm); xdm = fplmod(0, xmh, fsm);
                resh = [xdh(1); xdh(2); xdh(3); xdh(6); sqrt(xh(1)^2+xh(2)^2)-v_target];
                resm = [xdm(1); xdm(2); xdm(3); xdm(6); sqrt(xmh(1)^2+xmh(2)^2)-v_target];
                J(:,k) = (resh - resm) / (2*xs);
            end
            xt = xt - 0.7 * (J \ res);
        end
        x_tr(ivar) = xt;
        
        if ~converged
            fprintf('  V=%.1f H=%.0f: FAILED TO TRIM\n', v_target, h_target);
            continue;
        end
        
        % B. LINEARIZATION
        % Augmented State: [u w q theta v_int]
        state_idx = [1 2 3 4];
        n_aug = 5;
        A = zeros(n_aug, n_aug);
        xs_lin = 1e-5;
        
        cl_system = @(x_vec) closed_loop_dynamics(x_vec, x_tr, fsm, Kp_v, Ki_v, Kd_v, Kq, v_target);
        
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
        
        % Practical Stability: Accept real parts if Time-to-Double > 60s
        stability_threshold = log(2)/60;
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
set(gca, 'Color', [0.8 0.8 0.8]); % Gray background for NaN
imagesc(v_vec, h_vec, damping_matrix, 'AlphaData', ~isnan(damping_matrix));
set(gca, 'YDir', 'normal');
colorbar; colormap(gca, parula);
title('Damping Ratio \zeta (Higher is Better)');
xlabel('Airspeed (m/s)'); ylabel('Altitude (m)');
clim([0 1]); 
xticks(v_vec); % Show all airspeeds on X-axis

% Add Acceptable Envelope Contour (e.g., Zeta >= 0.3)
hold on;
[~, h_cont] = contour(V_grid, H_grid, damping_matrix, [0.3 0.3], 'r', 'LineWidth', 2.5, 'DisplayName', 'Acceptable Envelope (\zeta \geq 0.3)');
% No legend for imagesc/plot, only for h_cont
legend(h_cont, 'Location', 'southoutside');

subplot(1,2,2);
set(gca, 'Color', [0.8 0.8 0.8]); % Gray for NaN (Trim Failure)
intensity_map = t_half_matrix;
imagesc(v_vec, h_vec, intensity_map, 'AlphaData', ~isnan(intensity_map));
set(gca, 'YDir', 'normal');

% Custom Divergent Colormap
n_colors = 128;
cmap_g = [linspace(0.0, 0.8, n_colors/2)', linspace(0.3, 0.8, n_colors/2)', linspace(0.0, 0.8, n_colors/2)']; % Deeper Forest Green
cmap_r = [linspace(0.8, 0.6, n_colors/2)', linspace(0.8, 0.0, n_colors/2)', linspace(0.8, 0.0, n_colors/2)']; % Dark Crimson Red
cmap = [cmap_g; cmap_r];
colormap(gca, cmap);

max_viz = 0.5; % Slower saturation (0.5 instead of 0.1) softens the overall look
clim([-max_viz max_viz]);

cb = colorbar;
% Custom labels for the colorbar to show seconds
ticks = [-0.5 -0.2 -0.05 0 0.05 0.2 0.5];
tick_labels = cell(size(ticks));
for k = 1:length(ticks)
    if ticks(k) < -0.0001
        tick_labels{k} = sprintf('T_{1/2}=%.1fs', 1/abs(ticks(k)));
    elseif ticks(k) > 0.0001
        tick_labels{k} = sprintf('T_{2}=%.1fs', 1/ticks(k));
    else
        tick_labels{k} = 'Inf (Neutral)';
    end
end
set(cb, 'Ticks', ticks, 'TickLabels', tick_labels);

% --- ADD BLACK STRIPED HATCHING FOR NaN ---
for ax_idx = 1:2
    subplot(1,2,ax_idx);
    hold on;
    [row, col] = find(isnan(damping_matrix));
    for k = 1:length(row)
        % Cell bounds
        v = v_vec(col(k)); h = h_vec(row(k));
        dv = 2.5; dh = 250; % Half-widths
        % Parallel diagonal stripes (Black Striped)
        % Using plot with HandleVisibility off to keep legend clean
        plot([v-dv v], [h h+dh], 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot([v v+dv], [h-dh h], 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

title('Stability (Fast = Deep Color, 100s+ = Faint)');
xlabel('Airspeed (m/s)'); ylabel('Altitude (m)');
xticks(v_vec); % Show all airspeeds on X-axis
sgtitle('Elevator PID Performance: Airspeed-Hold Mode');
saveas(gcf, 'control_performance_envelope.png');

fprintf('\nAnalysis complete. Heatmaps saved.\n');

% --- HELPER: CLOSED LOOP DYNAMICS ---
function xdot_aug = closed_loop_dynamics(x_vec, x_tr, fsm, Kp, Ki, Kd, Kq, v_target)
    % x_vec: [u w q theta v_int]
    u = x_vec(1);
    w = x_vec(2);
    q = x_vec(3);
    theta = x_vec(4);
    v_int = x_vec(5);
    
    de_trim = x_tr(8);
    dp_trim = x_tr(9);
    
    v_curr = sqrt(u^2 + w^2);
    v_err = v_target - v_curr;
    
    % Get physics xdot for v_dot estimation
    xf_tmp = [u; w; q; theta; 0; x_tr(6); 0; de_trim; dp_trim];
    xd_tmp = fplmod(0, xf_tmp, fsm);
    v_dot = (u*xd_tmp(1) + w*xd_tmp(2)) / v_curr;
    
    % PID Logic (Speed on Elevator) - Kd correctly subtracted for damping
    de = de_trim + (Kp * v_err + Ki * v_int - Kd * v_dot) + Kq * q;
    
    % Full state for fplmod
    xf = [u; w; q; theta; 0; x_tr(6); 0; de; dp_trim];
    xd_phys = fplmod(0, xf, fsm);
    
    % Augmented derivative [udot wdot qdot thetadot v_err]
    xdot_aug = [xd_phys(1); xd_phys(2); xd_phys(3); xd_phys(4); v_err];
end
