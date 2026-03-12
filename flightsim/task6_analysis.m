% task6_analysis.m
% Comprehensive Longitudinal Stability and Control Analysis for the Windex
% 
% This script addresses Task 6 requirements:
% 1. Eigenvalue tracing (Root Locus vs Airspeed)

clear; close all;
fsm = make_fsim();

% --- PARAMETERS ---
v_sweep = 26:1:80; 
x_np = 0.43;
cref = fsm.cref;
sm_targets = [0.03, 0.06, 0.09, 0.12, 0.15];
colors = {'r', 'g', 'k', 'm', 'b'}; % Aft, 6%, Nom, 12%, Forward
labels = {'3% SM (Aft)', '6% SM', '9% SM (Nom)', '12% SM', '15% SM (Fwd)'};

all_results = struct();

for s = 1:length(sm_targets)
    current_sm = sm_targets(s);
    xcg = x_np - current_sm * cref;
    fsm.cog(1) = xcg;
    
    fprintf('\n--- Analyzing Static Margin: %.0f%% (x_cg = %.3f m) ---\n', current_sm*100, xcg);
    speed_data = [];
    
    % Starting state guess: [u w q theta dist alt fuel de dp]
    % For 3% SM (aft CG), start with less elevator deflection as requested
    if current_sm == 0.03
        x = [v_sweep(1) 0.3 0.0 0.05 0 1000 0 -0.05 0.1]'; 
    else
        x = [v_sweep(1) 0.4 0.0 0.05 0 1000 0 0.0 0.1]'; 
    end
    
    prev_eg = []; % For tracking eigenvalues
    
    for i = 1:length(v_sweep)
        vset = v_sweep(i);
        x(1) = vset; 
        
        ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8];
        xtrim = x(ivar);
        for iter = 1:50 
            x(ivar) = xtrim;
            [xdot] = fplmod(0, x, fsm);
            xstep = 1e-7;
            J = zeros(5, 5);
            for j = 1:5
                xh = x; xh(ivar(j)) = xh(ivar(j)) + xstep;
                xmh = x; xmh(ivar(j)) = xmh(ivar(j)) - xstep;
                [xdoth] = fplmod(0, xh, fsm);
                [xdotmh] = fplmod(0, xmh, fsm);
                J(:,j) = (xdoth(ifun) - xdotmh(ifun)) / (2*xstep);
            end
            ftrim = xdot(ifun); ftrim(5) = ftrim(5) - vset;
            if norm(ftrim) < 1e-9, break; end
            xtrim = xtrim - 0.7 * (J \ ftrim);
        end
        x(ivar) = xtrim;
        if norm(ftrim) > 1e-3, continue; end
        
        % Linearization
        st_idx = [1 2 3 4]; A = zeros(4,4); xstep = 1e-5;
        for j = 1:4
            idx = st_idx(j);
            xh = x; xh(idx) = xh(idx) + xstep;
            xmh = x; xmh(idx) = xmh(idx) - xstep;
            [xdoth] = fplmod(0, xh, fsm);
            [xdotmh] = fplmod(0, xmh, fsm);
            A(:,j) = (xdoth(st_idx) - xdotmh(st_idx)) / (2*xstep);
        end
        
        eg = eig(A);
        
        % Robust Mode Identification using displacement from previous speed
        if isempty(prev_eg)
            [~, s_idx] = sort(abs(real(eg)), 'descend');
            eg_sorted = eg(s_idx);
        else
            % Match to previous eigenvalues to prevent swapping (spikes)
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
        % Short Period (indices 1, 2)
        mode_info.t_fast = mean(t_vals(1:2));
        % Phugoid (indices 3, 4)
        mode_info.t_slow = mean(t_vals(3:4));
        
        v_idx = length(speed_data) + 1;
        speed_data(v_idx).speed = vset;
        speed_data(v_idx).eig = eg_sorted;
        speed_data(v_idx).mode = mode_info;
        speed_data(v_idx).de = x(8) * 180 / pi; % Elevator deflection in degrees
    end
    all_results(s).speed_data = speed_data;
end

% --- PLOTTING ---
figure(1); clf; set(gcf, 'Position', [100 100 1000 500], 'Color', 'w');
figure(2); clf; set(gcf, 'Position', [100 100 700 700], 'Color', 'w');
figure(3); clf; set(gcf, 'Position', [100 100 700 500], 'Color', 'w');

for s = 1:length(sm_targets)
    sd = all_results(s).speed_data;
    vs = [sd.speed];
    eigs = [sd.eig];
    md = [sd.mode];
    t_fast = [md.t_fast];
    t_slow = [md.t_slow];
    valid = ~isnan(t_fast) & ~isnan(t_slow);
    
    % Figure 1: Root Locus
    figure(1);
    subplot(1,2,1); hold on; grid on;
    plot(real(eigs), imag(eigs), [colors{s} '.'], 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(real(eigs(:,1)), imag(eigs(:,1)), [colors{s} 'o'], 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(real(eigs(:,end)), imag(eigs(:,end)), [colors{s} 'x'], 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % Marker proxy for legend
    plot(nan, nan, colors{s}, 'DisplayName', labels{s});
    
    subplot(1,2,2); hold on; grid on;
    plot(real(eigs(:,valid)), imag(eigs(:,valid)), [colors{s} '.'], 'MarkerSize', 8, 'HandleVisibility', 'off');
    plot(real(eigs(:,1)), imag(eigs(:,1)), [colors{s} 'o'], 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(real(eigs(:,end)), imag(eigs(:,end)), [colors{s} 'x'], 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Figure 2: Times
    figure(2);
    subplot(2,1,1); hold on; grid on;
    plot(vs(valid), t_fast(valid), [colors{s} '-s'], 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', labels{s});
    
    subplot(2,1,2); hold on; grid on;
    plot(vs(valid), t_slow(valid), [colors{s} '-o'], 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', labels{s});

    % Figure 3: Trimmed Elevator Deflection
    if isfield(sd, 'de')
        de = [sd.de];
        figure(3); hold on; grid on;
        plot(vs(valid), de(valid), [colors{s} '-^'], 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', labels{s});
    end
end

% Formatting Figure 1
figure(1);
subplot(1,2,1); title('Longitudinal Root Locus Comparison');
xlabel('Real Axis (1/s)'); ylabel('Imaginary Axis (rad/s)');
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
% Airspeed indicator markers in legend
plot(nan, nan, 'go', 'DisplayName', sprintf('Start: %.0f m/s', v_sweep(1)));
plot(nan, nan, 'rx', 'DisplayName', sprintf('End: %.0f m/s', v_sweep(end)));
legend('Location', 'best', 'FontSize', 9);

subplot(1,2,2); title('Phugoid Region Zoom');
xlabel('Real Axis (1/s)'); ylabel('Imaginary Axis (rad/s)');
line([0 0], [-0.5 0.5], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
axis([-0.06 0.06 -0.5 0.5]);

% Formatting Figure 2
figure(2);
subplot(2,1,1); title('Short Period Mode: Time to Half ($T_{1/2}$)', 'Interpreter', 'latex');
ylabel('Time (s) [Neg=-T_{1/2}]'); legend('Location', 'best', 'FontSize', 8);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

subplot(2,1,2); title('Phugoid Mode: Time to Double ($T_2$) or Half ($T_{1/2}$)', 'Interpreter', 'latex');
xlabel('Airspeed (m/s)'); ylabel('Time (s) [Pos=T_2, Neg=-T_{1/2}]');
ylim([-250 100]);
legend('Location', 'best', 'FontSize', 8);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

sgtitle('Stability Mode Evolution vs Static Margin');

% Formatting Figure 3
figure(3);
title('Trimmed Elevator Deflection vs Airspeed');
xlabel('Airspeed (m/s)'); ylabel('Elevator Deflection (deg)');
legend('Location', 'best', 'FontSize', 9);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

fprintf('\nSaving refined plots...\n');
saveas(1, 'task6_root_locus.png');
saveas(2, 'task6_time_to_half.png');
saveas(3, 'task6_elevator_trim.png');
fprintf('\nAnalysis complete.\n');
