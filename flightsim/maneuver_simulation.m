% maneuver_simulation.m
% Simulates a Stall and a Looping maneuver for the Windex aircraft.

clear; close all;
fsm = make_fsim();

%% --- 1. MANEUVER: STALL ---
fprintf('--- Maneuver 1: Stall Simulation ---\n');
v_stall_start = 30.0; % m/s
x0_stall = [v_stall_start 0.5 0.0 0.05 0 1000 0 0.0 0.0]'; 
x_trim_stall = trim_aircraft(x0_stall, v_stall_start, fsm);

% Stall Control Schedule
% Goal: Pitch up to Alpha > 10 deg, then recover.
t_stall = 0:0.1:20;
de_stall = x_trim_stall(8) * ones(size(t_stall));
% Pull up at t=2s
de_stall(t_stall >= 2 & t_stall < 10) = x_trim_stall(8) - 5*pi/180; % Large nose-up
% Release/Recover at t=7s
% de_stall(t_stall >= 5 & t_stall < 7) = x_trim_stall(8) + 5*pi/180; % Nose-down recovery
% Return to trim
de_stall(t_stall >= 10) = x_trim_stall(8);

fsm.de_set = [t_stall' de_stall'];
fsm.dp_set = x_trim_stall(9);

fprintf('Running Stall Simulation...\n');
[t_vec_s, Y_s] = ode45(@(t,y) simsub(t, y, fsm), [0 20], x_trim_stall(1:7));

% Calculate Alpha and Nz for Stall
alpha_s = zeros(length(t_vec_s), 1);
nz_s = zeros(length(t_vec_s), 1);
for i = 1:length(t_vec_s)
    [xdot] = fplmod(t_vec_s(i), [Y_s(i,:) interp1(t_stall, de_stall, t_vec_s(i)) fsm.dp_set]', fsm);
    alpha_s(i) = xdot(9);
    nz_s(i) = xdot(13);
end

%% --- 2. MANEUVER: LOOPING ---
fprintf('\n--- Maneuver 2: Looping Simulation ---\n');
v_loop_start = 75.0; % m/s
x0_loop = [v_loop_start 0.5 0.0 0.05 0 1000 0 0.0 0.1]'; % Start with some thrust if needed
x_trim_loop = trim_aircraft(x0_loop, v_loop_start, fsm);

% Looping Control Schedule
% Goal: Pull into a loop without exceeding 6g.
t_loop = 0:0.1:15;
de_loop = x_trim_loop(8) * ones(size(t_loop));
% Pull up hard at t=1s
de_loop(t_loop >= 1 & t_loop < 12.5) = x_trim_loop(8) - 3.5*pi/180; 

de_loop(t_loop >= 12.5) = x_trim_loop(8);


fsm.de_set = [t_loop' de_loop'];
fsm.dp_set = x_trim_loop(9);

fprintf('Running Looping Simulation...\n');
[t_vec_l, Y_l] = ode45(@(t,y) simsub(t, y, fsm), [0 15], x_trim_loop(1:7));

% Calculate Alpha and Nz for Loop
alpha_l = zeros(length(t_vec_l), 1);
nz_l = zeros(length(t_vec_l), 1);
for i = 1:length(t_vec_l)
    [xdot] = fplmod(t_vec_l(i), [Y_l(i,:) interp1(t_loop, de_loop, t_vec_l(i)) fsm.dp_set]', fsm);
    alpha_l(i) = xdot(9);
    nz_l(i) = xdot(13);
end

%% --- PLOTTING ---
figure(1); clf; set(gcf, 'Color', 'w', 'Name', 'Stall Maneuver');
subplot(3,1,1); hold on;
plot(t_vec_s, alpha_s*180/pi, 'b', 'LineWidth', 1.5, 'DisplayName', 'Alpha');
plot(t_vec_s, interp1(t_stall, de_stall, t_vec_s)*180/pi, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Elevator');
grid on; ylabel('Angle (deg)'); title('Maneuver 1: Stall (V_{start}=30m/s)');
line(xlim, [10 10], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off'); % Stall target
legend('Location', 'best');

subplot(3,1,2); plot(t_vec_s, Y_s(:,6), 'k', 'LineWidth', 1.5); grid on;
ylabel('Altitude (m)');

subplot(3,1,3); hold on;
plot(t_vec_s, Y_s(:,1), 'r', 'LineWidth', 1.2, 'DisplayName', 'u (Forward)');
plot(t_vec_s, Y_s(:,2), 'g', 'LineWidth', 1.2, 'DisplayName', 'w (Downward)');
plot(t_vec_s, sqrt(Y_s(:,1).^2 + Y_s(:,2).^2), 'b--', 'LineWidth', 1.2, 'DisplayName', 'V Total');
grid on; ylabel('Velocity (m/s)'); xlabel('Time (s)');
legend('Location', 'best'); title('Velocity Components during Stall');

figure(2); clf; set(gcf, 'Color', 'w', 'Name', 'Looping Maneuver');
subplot(2,2,1); plot(Y_l(:,5), Y_l(:,6), 'b', 'LineWidth', 1.5); grid on;
xlabel('Distance (m)'); ylabel('Altitude (m)'); title('Looping Trajectory');
axis equal;

subplot(2,2,2); plot(t_vec_l, nz_l, 'r', 'LineWidth', 1.5); grid on;
ylabel('G-load (Total Magnitude)'); title('Total G-force during Loop');
line(xlim, [6 6], 'Color', 'k', 'LineStyle', '--'); % 6g Limit

subplot(2,2,3); plot(t_vec_l, Y_l(:,1), 'g', 'LineWidth', 1.5); grid on;
ylabel('Velocity u (m/s)'); xlabel('Time (s)');

subplot(2,2,4); hold on;
plot(t_vec_l, alpha_l*180/pi, 'm', 'LineWidth', 1.5, 'DisplayName', 'Alpha');
plot(t_vec_l, interp1(t_loop, de_loop, t_vec_l)*180/pi, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Elevator');
grid on; ylabel('Angle (deg)'); xlabel('Time (s)');
legend('Location', 'best');

fprintf('\nSimulations Complete. Review the plots.\n');

%% --- HELPER: TRIM FUNCTION ---
function x = trim_aircraft(x0, vset, fsm)
    ivar = [1 2 4 8 9]; % u, w, theta, delta_e, delta_p
    ifun = [1 2 3 6 8]; % udot, wdot, qdot, altdot (for glide), airspeed
    
    x = x0;
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
        
        ftrim = xdot(ifun);
        ftrim(5) = ftrim(5) - vset; 
        
        if norm(ftrim) < 1e-9, break; end
        xtrim = xtrim - 0.7 * (J \ ftrim);
    end
    x(ivar) = xtrim;
end
