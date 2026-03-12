% autothrottle_design.m
% PID Autothrottle with Altitude Hold and Coordinated Feedforward
clear; close all; clc;
fsm = make_fsim();

% --- 1. ROBUST TRIM SWEEP (30 m/s -> [24, 36]) ---
% We sweep outwards from 30 m/s to ensure convergence
v_points = [30 29 28 27 26 25 24 31 32 33 34 35 36];
v_table = sort(v_points);
de_table = zeros(size(v_table));
dp_table = zeros(size(v_table));

fprintf('Calculating robust trim table for speed scheduling...\n');
x_guess = [30 0.5 0 0.05 0 1000 0 0.02 0.08]'; % Good initial guess for 30
for v_test = [30:-1:24, 31:1:36]
    xt = x_guess([1 2 4 8 9]);
    for iter = 1:50
        x_guess([1 2 4 8 9]) = xt;
        xd = fplmod(0, x_guess, fsm);
        xs = 1e-8; J = zeros(5, 5);
        for j = 1:5
            xh = x_guess; vars=[1 2 4 8 9]; xh(vars(j)) = xh(vars(j)) + xs;
            xmh = x_guess; xmh(vars(j)) = xmh(vars(j)) - xs;
            xd_h = fplmod(0, xh, fsm); xd_m = fplmod(0, xmh, fsm);
            J(:,j) = (xd_h([1 2 3 6 8]) - xd_m([1 2 3 6 8])) / (2*xs);
        end
        res = xd([1 2 3 6 8]); res(5) = res(5) - v_test;
        if norm(res) < 1e-11, break; end
        xt = xt - 0.7 * (J \ res);
    end
    idx = find(v_table == v_test);
    de_table(idx) = x_guess(8);
    dp_table(idx) = x_guess(9);
    % x_guess stays updated for next point in sweep
end

% Reference values for 30 m/s
de_trim_30 = interp1(v_table, de_table, 30);
dp_trim_30 = interp1(v_table, dp_table, 30);

fprintf('\n--- System Parameters ---\n');
fprintf('100%% Throttle = %.1f N Thrust\n', fsm.thrustmax);
fprintf('Trim at 30 m/s: Throttle = %.2f%%, Elevator = %.2f deg\n', dp_trim_30*100, de_trim_30*180/pi);
fprintf('Limited Schedule [V : de]:\n');
for i=find(v_table>=26 & v_table<=34)
    fprintf('  %2.0f m/s : %6.2f deg\n', v_table(i), de_table(i)*180/pi);
end

% --- 2. CONTROL GAINS ---
Kp = 0.015; Ki = 0.002;  % Speed -> Throttle
Kh = 0.005; Kq = 0.5;   % Altitude/Pitch -> Elevator (Stabilization)

% --- 3. SIMULATION ---
dt = 0.05; t_end = 150; t = 0:dt:t_end; t_nudge = 10;
scenarios = {'NudgeDown29', 'NudgeUp31'};
results = struct();

for s = 1:2
    mode = scenarios{s};
    x = zeros(8, 1); % [u w q theta x h fuel x_int]
    
    % Initialize precisely at 30 m/s trim
    x_tr = [30 0.5 0 0.05 0 1000 0 3.85*pi/180 0.0827]';
    xt = x_tr([1 2 4 8 9]);
    for it=1:50
        x_tr([1 2 4 8 9]) = xt; xd = fplmod(0,x_tr,fsm);
        res = xd([1 2 3 6 8]); res(5) = res(5)-30;
        if norm(res) < 1e-12, break; end
        xs=1e-8; J=zeros(5,5); 
        for j=1:5
            xh=x_tr; vars=[1 2 4 8 9]; xh(vars(j))=xh(vars(j))+xs; 
            xmh=x_tr; vars=[1 2 4 8 9]; xmh(vars(j))=xmh(vars(j))-xs;
            xd_h = fplmod(0,xh,fsm); xd_m = fplmod(0,xmh,fsm);
            J(:,j)=(xd_h([1 2 3 6 8])-xd_m([1 2 3 6 8]))/(2*xs);
        end
        xt = xt - 0.7*(J\res);
    end
    x(1:7) = x_tr(1:7); x(8) = 0;
    
    nudge_done = false;
    history_u = zeros(length(t), 1); history_dp = zeros(length(t), 1);
    history_de = zeros(length(t), 1); history_h = zeros(length(t), 1);
    
    for i = 1:length(t)
        curr_t = t(i);
        curr_u = sqrt(x(1)^2 + x(2)^2);
        
        if curr_t >= t_nudge && ~nudge_done
            if strcmp(mode, 'NudgeDown29'), x(1) = 29.0; else, x(1) = 31.0; end
            nudge_done = true; x(8) = 0; curr_u = sqrt(x(1)^2 + x(2)^2);
        end
        
        [~, de, dp] = control_wrapper(curr_t, x, fsm, v_table, de_table, dp_trim_30, Kp, Ki, Kh, Kq, 30.0, 1000.0);
        history_u(i) = curr_u; history_dp(i) = dp; history_de(i) = de; history_h(i) = x(6);
        
        k1 = control_wrapper(curr_t, x, fsm, v_table, de_table, dp_trim_30, Kp, Ki, Kh, Kq, 30.0, 1000.0);
        k2 = control_wrapper(curr_t + dt/2, x + k1*dt/2, fsm, v_table, de_table, dp_trim_30, Kp, Ki, Kh, Kq, 30.0, 1000.0);
        k3 = control_wrapper(curr_t + dt/2, x + k2*dt/2, fsm, v_table, de_table, dp_trim_30, Kp, Ki, Kh, Kq, 30.0, 1000.0);
        k4 = control_wrapper(curr_t + dt, x + k3*dt, fsm, v_table, de_table, dp_trim_30, Kp, Ki, Kh, Kq, 30.0, 1000.0);
        x = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    results.(mode).u = history_u; results.(mode).dp = history_dp;
    results.(mode).de = history_de; results.(mode).h = history_h;
end

% --- 4. VISUALIZATION (With Axis Limits) ---
figure(1); clf; set(gcf, 'Position', [100 100 1100 850], 'Color', 'w');
subplot(3,2,1); hold on; grid on; title('Speed Recovery (Nudge to 29)');
plot(t, results.NudgeDown29.u, 'b', 'LineWidth', 1.5); ylim([25 35]); ylabel('u (m/s)');
subplot(3,2,3); hold on; grid on; title('Throttle Output');
plot(t, results.NudgeDown29.dp, 'r', 'LineWidth', 1.2); ylim([0 1]); ylabel('dp (0-1)');
subplot(3,2,5); hold on; grid on; title('Elevator (Trim + Alt-Hold + Damp)');
plot(t, results.NudgeDown29.de*180/pi, 'g', 'LineWidth', 1.2); ylim([-5 15]); ylabel('de (deg)'); xlabel('Time (s)');

subplot(3,2,2); hold on; grid on; title('Speed Recovery (Nudge to 31)');
plot(t, results.NudgeUp31.u, 'b', 'LineWidth', 1.5); ylim([25 45]); ylabel('u (m/s)');
subplot(3,2,4); hold on; grid on; title('Throttle Output');
plot(t, results.NudgeUp31.dp, 'r', 'LineWidth', 1.2); ylim([0 1]); ylabel('dp (0-1)');
subplot(3,2,6); hold on; grid on; title('Elevator (Trim + Alt-Hold + Damp)');
plot(t, results.NudgeUp31.de*180/pi, 'g', 'LineWidth', 1.2); ylim([-5 15]); ylabel('de (deg)'); xlabel('Time (s)');
sgtitle('Autothrottle Coordinated Control - Windex Aircraft');
saveas(1, 'autothrottle_coordinated_results.png');

figure(2); clf; set(gcf, 'Position', [100 100 800 600], 'Color', 'w');
subplot(2,1,1); hold on; grid on; title('Altitude Stability (Nudge Down)');
plot(t, results.NudgeDown29.h - 1000, 'm', 'LineWidth', 1.5); ylim([-20 20]); ylabel('\Delta h (m)');
subplot(2,1,2); hold on; grid on; title('Altitude Stability (Nudge Up)');
plot(t, results.NudgeUp31.h - 1000, 'm', 'LineWidth', 1.5); ylim([-20 20]); ylabel('\Delta h (m)'); xlabel('Time (s)');
sgtitle('Altitude Maintenance via Scheduled Trim & Active Feedback');
saveas(2, 'autothrottle_coordinated_coupling.png'); fprintf('\nSimulations complete.\n');

function [xdot, de, dp] = control_wrapper(t, x, fsm, v_table, de_table, dp_trim, Kp, Ki, Kh, Kq, v_target, h_target)
    spd = sqrt(x(1)^2 + x(2)^2); h_err = x(6) - h_target;
    u_err = v_target - spd;
    dp = max(0, min(1, dp_trim + Kp*u_err + Ki*x(8)));
    de = interp1(v_table, de_table, spd, 'linear', 'extrap') + Kh*h_err + Kq*x(3);
    xf = zeros(13, 1); xf(1:7) = x(1:7); xf(8) = de; xf(9) = dp;
    xdot_phys = fplmod(t, xf, fsm);
    xdot = zeros(8, 1); xdot(1:7) = xdot_phys(1:7); xdot(8) = u_err;
end
