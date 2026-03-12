% task6_glide_analysis.m
% Longitudinal Stability Analysis for an UNPOWERED Windex (Engine Out)
% 
% This script traces Eigenvalues vs Airspeed, exactly like task6_analysis.m,
% but it forces the throttle to zero (dp = 0) and allows the aircraft to
% descend (altdot ~= 0) to maintain trim, completely eliminating the 
% pylon engine's nose-down pitching moment.

clear; close all;
addpath('flightsim');
fsm = make_fsim();

% --- PARAMETERS ---
v_sweep = 26:1:80; 
x_np = 0.43;
cref = fsm.cref;
sm_targets = [0.03, 0.06, 0.09, 0.12, 0.15];
colors = {'r', 'g', 'k', 'm', 'b'};
labels = {'3% SM (Aft)', '6% SM', '9% SM (Nom)', '12% SM', '15% SM (Fwd)'};

all_results = struct();

fprintf('=== Unpowered Glide Stability Analysis (Engine Out) ===\n');

for s = 1:length(sm_targets)
    current_sm = sm_targets(s);
    xcg = x_np - current_sm * cref;
    fsm.cog(1) = xcg;
    
    fprintf('\n--- Analyzing Glide at Static Margin: %.0f%% ---\n', current_sm*100);
    speed_data = [];
    
    % Starting guess
    % Notice we only guess 4 variables to change (u, w, theta, de)
    x = [v_sweep(1) 0.8 0.0 -0.05 0 1000 0 -0.0 0.0]'; 
    
    prev_eg = [];
    
    for i = 1:length(v_sweep)
        vset = v_sweep(i);
        x(1) = vset; 
        x(9) = 0; % FORCED ZERO THROTTLE (Engine Out)
        
        % GLIDE TRIM CONFIGURATION
        % We do NOT solve for dp, and we do NOT enforce altdot = 0
        ivar = [1 2 4 8]; % u, w, theta, de
        ifun = [1 2 3 8]; % udot, wdot, qdot, airspeed
        
        xtrim = x(ivar);
        for iter = 1:50 
            x(ivar) = xtrim;
            [xdot] = fplmod(0, x, fsm);
            
            xstep = 1e-6;
            J = zeros(4, 4);
            for j = 1:4
                xh = x; xh(ivar(j)) = xh(ivar(j)) + xstep;
                xmh = x; xmh(ivar(j)) = xmh(ivar(j)) - xstep;
                [xdoth] = fplmod(0, xh, fsm);
                [xdotmh] = fplmod(0, xmh, fsm);
                J(:,j) = (xdoth(ifun) - xdotmh(ifun)) / (2*xstep);
            end
            
            ftrim = xdot(ifun); 
            ftrim(4) = ftrim(4) - vset; % Enforce target airspeed
            
            if norm(ftrim) < 1e-7, break; end
            xtrim = xtrim - 0.7 * (J \ ftrim);
        end
        x(ivar) = xtrim;
        if norm(ftrim) > 1e-3
            fprintf('Warning: Failed to trim at V = %.1f m/s\n', vset);
            continue; 
        end
        
        % Linearization arrays (4 longitudinal states: u, w, q, theta)
        st_idx = [1 2 3 4]; 
        A = zeros(4,4); 
        xstep = 1e-5;
        for j = 1:4
            idx = st_idx(j);
            xh = x; xh(idx) = xh(idx) + xstep;
            xmh = x; xmh(idx) = xmh(idx) - xstep;
            [xdoth] = fplmod(0, xh, fsm);
            [xdotmh] = fplmod(0, xmh, fsm);
            A(:,j) = (xdoth(st_idx) - xdotmh(st_idx)) / (2*xstep);
        end
        
        eg = eig(A);
        
        % Robust Mode Identification
        if isempty(prev_eg)
            [~, s_idx] = sort(abs(real(eg)), 'descend');
            eg_sorted = eg(s_idx);
        else
            new_idx = zeros(4,1);
            available = 1:4;
            for k = 1:4
                dist = abs(eg(available) - prev_eg(k));
                [~, best_loc] = min(dist);
                best_idx = available(best_loc);
                new_idx(k) = best_idx;
                available(best_loc) = [];
            end
            eg_sorted = eg(new_idx);
        end
        prev_eg = eg_sorted;
        
        % Time to half/double values
        t_vals = log(2) ./ real(eg_sorted);
        
        mode_info = struct();
        mode_info.t_fast = mean(t_vals(1:2)); % Short Period
        mode_info.t_slow = mean(t_vals(3:4)); % Phugoid
        
        v_idx = length(speed_data) + 1;
        speed_data(v_idx).speed = vset;
        speed_data(v_idx).eig = eg_sorted;
        speed_data(v_idx).mode = mode_info;
    end
    all_results(s).speed_data = speed_data;
end

% --- PLOTTING ---
figure(1); clf; set(gcf, 'Position', [100 100 700 700], 'Color', 'w');

for s = 1:length(sm_targets)
    sd = all_results(s).speed_data;
    if isempty(sd), continue; end
    vs = [sd.speed];
    md = [sd.mode];
    t_fast = [md.t_fast];
    t_slow = [md.t_slow];
    valid = ~isnan(t_fast) & ~isnan(t_slow);
    
    subplot(2,1,1); hold on; grid on;
    plot(vs(valid), t_fast(valid), [colors{s} '-s'], 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', labels{s});
    
    subplot(2,1,2); hold on; grid on;
    plot(vs(valid), t_slow(valid), [colors{s} '-o'], 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', labels{s});
end

figure(1);
subplot(2,1,1); title('Short Period Mode: Unpowered Glide ($T_{1/2}$)', 'Interpreter', 'latex');
ylabel('Time (s) [Neg=-T_{1/2}]'); legend('Location', 'best', 'FontSize', 8);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

subplot(2,1,2); title('Phugoid Mode: Unpowered Glide ($T_2$ or $T_{1/2}$)', 'Interpreter', 'latex');
xlabel('Airspeed (m/s)'); ylabel('Time (s) [Pos=T_2, Neg=-T_{1/2}]');
ylim([-250 100]);
legend('Location', 'best', 'FontSize', 8);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

sgtitle('Unpowered Stability Mode Evolution vs Static Margin (Engine Out)');

fprintf('\nSaving plots...\n');
saveas(1, 'task6_glide_stability.png');
fprintf('Done.\n');
