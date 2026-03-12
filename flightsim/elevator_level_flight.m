% elevator_level_flight.m
% Elevator PID Controller for Level Flight (Altitude Hold)
% This script demonstrates how the elevator compensates for throttle changes
% to maintain a constant altitude.

clear; close all; clc;
fsm = make_fsim();

% --- 1. INITIAL TRIM at 30 m/s ---
fprintf('Trimming aircraft for 30 m/s level flight...\n');
v_start = 30.0;
h_target = 1000.0;
% [u w q theta x h fuel de dp]
x_trim = [v_start 0.5 0 0.05 0 h_target 0 0.06 0.08]';

% Simple robust trim solver
ivar = [1 2 4 8 9]; % u, w, theta, de, dp
ifun = [1 2 3 6 8]; % udot, wdot, qdot, hdot, airspeed
xt = x_trim(ivar);
for iter = 1:50
    x_trim(ivar) = xt;
    xd = fplmod(0, x_trim, fsm);
    res = xd(ifun); res(5) = res(5) - v_start;
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
% Altitude PID on Elevator
Kp_h = 0.005;   % Altitude error to elevator rad/m
Ki_h = 0.0002;  % Integral gain
Kd_h = 0.015;   % Altitude rate damping rad/(m/s)
Kq   = 0.4;     % Pitch rate damping rad/(rad/s)

% --- 3. SIMULATION ---
dt = 0.05;
t_end = 100;
t = 0:dt:t_end;

% Initialize state
x = x_trim(1:7);
h_int = 0;

% Storage
history = struct();
history.u = zeros(size(t));
history.h = zeros(size(t));
history.de = zeros(size(t));
history.dp = zeros(size(t));
history.theta = zeros(size(t));

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
    
    % PID CONTROL FOR ELEVATOR
    h_err = h_target - x(6);
    h_int = h_int + h_err * dt;
    h_dot = - ( -sin(x(4))*x(1) + cos(x(4))*x(2)); % Current vertical speed (up is positive)
    
    % The classic PID for elevator (Altitude Hold)
    % de_cmd = de_trim + feedback
    % Note: h_err positive (too low) -> want nose up -> decrease de (more negative/less positive?)
    % Wait, fplmod convention: +de is trailing edge down (nose down).
    % So too low (h_err > 0) -> need nose up -> -de.
    de = de_trim - (Kp_h * h_err + Ki_h * h_int - Kd_h * h_dot) + Kq * x(3);
    
    % saturation
    de = max(-20*pi/180, min(10*pi/180, de));
    
    % Record history
    history.u(i) = sqrt(x(1)^2 + x(2)^2);
    history.h(i) = x(6);
    history.de(i) = de;
    history.dp(i) = dp;
    history.theta(i) = x(4);
    
    % Integration (RK4)
    k1 = dxdt(x, de, dp, fsm);
    k2 = dxdt(x + k1*dt/2, de, dp, fsm);
    k3 = dxdt(x + k2*dt/2, de, dp, fsm);
    k4 = dxdt(x + k3*dt, de, dp, fsm);
    x = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% --- 4. PLOTTING ---
figure('Name', 'Elevator PID Level Flight', 'Color', 'w', 'Position', [100 100 1000 800]);

subplot(4,1,1);
plot(t, history.dp * 100, 'r', 'LineWidth', 1.5);
grid on; ylabel('Throttle (%)'); title('User Throttle Input');

subplot(4,1,2);
plot(t, history.h, 'b', 'LineWidth', 1.5);
hold on; line(xlim, [h_target h_target], 'Color', 'k', 'LineStyle', '--');
grid on; ylabel('Altitude (m)'); title('Altitude Maintenance (Level Flight)');

subplot(4,1,3);
plot(t, history.de * 180/pi, 'g', 'LineWidth', 1.5);
grid on; ylabel('Elevator (deg)'); title('PID Elevator Response');

subplot(4,1,4);
yyaxis left
plot(t, history.u, 'k', 'LineWidth', 1.5);
ylabel('Airspeed (m/s)');
yyaxis right
plot(t, history.theta * 180/pi, 'm', 'LineWidth', 1.2);
ylabel('Pitch (deg)');
grid on; xlabel('Time (s)');
legend('Airspeed', 'Pitch');

sgtitle('Windex Elevator PID: Maintaining Level Flight during Throttle Changes');
saveas(gcf, 'elevator_pid_results.png');

fprintf('Simulation complete. Plot saved as elevator_pid_results.png\n');

function xdot = dxdt(x, de, dp, fsm)
    xv = [x; de; dp];
    res = fplmod(0, xv, fsm);
    xdot = res(1:7);
end
