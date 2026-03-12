% Elevator PID Controller for Airspeed Hold
% This script demonstrates how the elevator maintains a constant airspeed
% in response to throttle changes by adjusting pitch.

clear; close all; clc;
fsm = make_fsim();

fprintf('Trimming aircraft for 30 m/s level flight...\n');
v_target = 30.0;
h_start = 1000.0;
% [u w q theta x h fuel de dp]
x_trim = [v_target 0.5 0 0.05 0 h_start 0 0.06 0.08]';

% Simple robust trim solver
ivar = [1 2 4 8 9]; % u, w, theta, de, dp
ifun = [1 2 3 6 8]; % udot, wdot, qdot, hdot, airspeed
xt = x_trim(ivar);
for iter = 1:50
    x_trim(ivar) = xt;
    xd = fplmod(0, x_trim, fsm);
    res = xd(ifun); res(5) = res(5) - v_target;
    if norm(res) < 1e-11, break; end
    xs = 1e-8; J = zeros(5,5);
    for j = 1:5
        xh = x_trim; xh(ivar(j)) = xh(ivar(j)) + xs;
        xmh = x_trim; xmh(ivar(j)) = xmh(ivar(j)) - xs;
        xdh = fplmod(0, xh, fsm); xdm = fplmod(0, xmh, fsm);
        J(:,j) = (xdh(ifun) - xdm(ifun)) / (2*xs);
    end
    xt = xt - 0.7 * (J \ res);
end
x_trim(ivar) = xt;
de_trim = x_trim(8);
dp_trim = x_trim(9);

fprintf('Trim complete: Throttle = %.2f%%, Elevator = %.2f deg\n', dp_trim*100, de_trim*180/pi);

% --- 2. CONTROL CONSTANTS ---
% Airspeed PID on Elevator (Optimized for Smoothed Cruise)
Kp_v = 0.001785;
Ki_v = 0.000259;
Kd_v = 0.018147;
% Actuator Constraints
de_rate_limit = 60 * pi/180; % 60 deg/sec max rate
tau_act = 0.1;           % Actuator time constant (s) - adds realistic lag
Kq = 0.054794; % Optimized smoothed damping
de_actual = de_trim; 
de_prev = de_trim;   % For rate limiting logic (initialized after trim)


% --- 3. SIMULATION ---
dt = 0.05;
t_end = 100;
t = 0:dt:t_end;

% Initialize state
x = x_trim(1:7);
v_int = 0;
% Initialize Storage History
history = struct();
history.u = zeros(size(t));
history.h = zeros(size(t));
history.de = zeros(size(t));
history.de_cmd = zeros(size(t));
history.dp = zeros(size(t));
history.alpha = zeros(size(t));

fprintf('Starting simulation with throttle steps...\n');

for i = 1:length(t)
    curr_t = t(i);
    
    % USER THROTTLE INPUT (The "User wants" part)
    % Initial steady for 5s, then ramp up throttle
    if curr_t < 10
        dp = dp_trim;
    elseif curr_t < 40
        % Linear ramp up to +50% throttle over 30s
        dp = dp_trim + 0.1 * (curr_t - 10)/30;
    elseif curr_t < 70
        % Step back down
        dp = dp_trim + 0.05;
    else
        % Return to trim throttle
        dp = dp_trim;
    end
    
    % PID CONTROL FOR ELEVATOR (Airspeed Hold)
    v_curr = sqrt(x(1)^2 + x(2)^2);
    v_err = v_target - v_curr;
    v_int = v_int + v_err * dt;
    
    % Estimate v_dot using udot and wdot (simplified)
    xd_tmp = fplmod(0, [x; 0; dp], fsm);
    v_dot = (x(1)*xd_tmp(1) + x(2)*xd_tmp(2)) / v_curr;
    
    % de_cmd = de_trim + feedback
    % v_err > 0 (Too slow) -> Need nose down -> +de
    % --- Actuator Dynamics (Lag + Rate Limit) ---
    % 1. Command from PID
    de_cmd = de_trim + (Kp_v * v_err + Ki_v * v_int - Kd_v * v_dot) + Kq * x(3);
    
    % 2. First-order lag (Servo response)
    % de_dot = (de_cmd - de_actual) / tau_act
    de_dot_cmd = (de_cmd - de_actual) / tau_act;
    
    % 3. Apply physical Rate Limit
    de_dot = max(-de_rate_limit, min(de_rate_limit, de_dot_cmd));
    
    % 4. Update actual position
    de_actual = de_actual + de_dot * dt;
    de = de_actual;
    
    % record the command too for comparison
    history.de_cmd(i) = de_cmd;
    
    % saturation
    de = max(-20*pi/180, min(10*pi/180, de));
    
    % Record history
    history.u(i) = sqrt(x(1)^2 + x(2)^2);
    history.h(i) = x(6);
    history.de(i) = de;
    history.dp(i) = dp;
    history.alpha(i) = atan2(x(2), x(1));
    
    % Integration (RK4)
    k1 = dxdt(x, de, dp, fsm);
    k2 = dxdt(x + k1*dt/2, de, dp, fsm);
    k3 = dxdt(x + k2*dt/2, de, dp, fsm);
    k4 = dxdt(x + k3*dt, de, dp, fsm);
    x = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% --- 4. PLOTTING ---
figure('Name', 'Elevator PID Level Flight', 'Color', 'w', 'Position', [100 100 1000 800]);

% 1. Throttle Input
subplot(3,1,1);
plot(t, history.dp * 100, 'r', 'LineWidth', 1.5);
grid on; ylabel('Throttle (%)'); title('User Throttle Input');
set(gca, 'YColor', 'r');

% 2. Airspeed (L) and Altitude (R)
subplot(3,1,2);
yyaxis left
plot(t, history.u, 'b', 'LineWidth', 1.5); hold on;
yline(v_target, 'b--', 'LineWidth', 1, 'HandleVisibility', 'off');
ylabel('Airspeed (m/s)'); set(gca, 'YColor', 'b');
yyaxis right
plot(t, history.h, 'k', 'LineWidth', 1.5);
ylabel('Altitude (m)'); set(gca, 'YColor', 'k');
grid on; title('Speed Maintenance and Altitude Response');

% 3. Elevator (L) and Angle of Attack (R)
subplot(3,1,3);
yyaxis left
plot(t, history.de_cmd * 180/pi, 'r--', 'LineWidth', 1); hold on;
plot(t, history.de * 180/pi, 'b', 'LineWidth', 1.5);
ylabel('Elevator (deg)'); set(gca, 'YColor', 'b');
yyaxis right
plot(t, history.alpha * 180/pi, 'm', 'LineWidth', 1.2);
ylabel('Angle of Attack (deg)'); set(gca, 'YColor', 'm');
grid on; title('Controller Effort and Angle of Attack');
xlabel('Time (s)');
legend('Elev Cmd', 'Elev Actual', 'AoA', 'Location', 'best');

sgtitle('Windex Elevator PID: Compact Control Performance');
saveas(gcf, 'elevator_pid_results.png');

fprintf('Simulation complete. Plot saved as elevator_pid_results.png\n');

function xdot = dxdt(x, de, dp, fsm)
    xv = [x; de; dp];
    res = fplmod(0, xv, fsm);
    xdot = res(1:7);
end
