% explore_controls.m
% 3D Control Exploration Script for the Windex.
% Analyzes the time-domain response to Elevator, Aileron, Rudder, and Flap inputs.

clear; close all;
fsm = make_fsim();

% 1. Trim at a target speed in 3D (GLIDER MODE)
V_target = 30.0;
x3d = zeros(16,1);
% Initial guess: [u v w p q r phi theta psi xe ye alt de da dr dp]
x3d(1) = V_target;
x3d(12) = 500;   % Start at 500m
x3d(13) = -0.05;  % Elevator guess
x3d(16) = 0.0;    % ZERO THRUST (Glider)

% 4-Variable Glide Trim Loop: [u w theta de]
% We solve for equilibrium at zero thrust. The plane will sink (altdot < 0).
iv = [1 3 8 13]; 
fprintf('Trimming glider for V = %.1f m/s (Zero Thrust)...\n', V_target);

for iter = 1:15
    [xdot] = fplmod3d(0, x3d, fsm);
    
    % Physics Targets: udot=0, wdot=0, qdot=0
    f = [xdot(1); xdot(3); xdot(5)];
    % Speed Target: V = V_target
    V_now = sqrt(x3d(1)^2 + x3d(2)^2 + x3d(3)^2);
    f(4) = V_now - V_target;
    
    if norm(f) < 1e-8
        fprintf('Glide trim converged in %d iterations.\n', iter);
        break;
    end
    
    % Jacobian Calculation
    J = zeros(4,4); xstep = 1e-6;
    for j = 1:4
        xh = x3d; xh(iv(j)) = xh(iv(j)) + xstep;
        xmh = x3d; xmh(iv(j)) = xmh(iv(j)) - xstep;
        
        [xdoth] = fplmod3d(0, xh, fsm);
        [xdotmh] = fplmod3d(0, xmh, fsm);
        
        fh = [xdoth(1); xdoth(3); xdoth(5)];
        fh(4) = sqrt(xh(1)^2+xh(2)^2+xh(3)^2) - V_target;
        
        fm = [xdotmh(1); xdotmh(3); xdotmh(5)];
        fm(4) = sqrt(xmh(1)^2+xmh(2)^2+xmh(3)^2) - V_target;
        
        J(:,j) = (fh - fm) / (2*xstep);
    end
    
    % Update guess
    x3d(iv) = x3d(iv) - (J \ f);
end

% --- CHOOSE A SCENARIO ---
scenario = 'aileron_roll'; % Options: 'elevator_pulse', 'aileron_roll', 'rudder_kick', 'flap_deployment'

t_end = 60; 
fsm.ctrl_set = zeros(6, 5); % [Time, de, da, dr, dp]
fsm.ctrl_set(:,1) = [0 2 3 5 6 t_end]';
fsm.ctrl_set(:,2) = x3d(13); % de (Trim)
fsm.ctrl_set(:,3) = 0;       % da (None)
fsm.ctrl_set(:,4) = 0;       % dr (None)
fsm.ctrl_set(:,5) = x3d(16); % dp (Trim)
fsm.delta_flap_static = 0;   % Initial flaps (rad)

switch scenario
    case 'elevator_pulse'
        mag = -2.0 * pi/180; % 2 deg up
        fsm.ctrl_set(3:4, 2) = x3d(13) + mag;
        t_end = 180; fsm.ctrl_set(end,1) = t_end; % Long time for Phugoid
        
    case 'aileron_roll'
        mag = 10.0 * pi/180; % 10 deg right
        fsm.ctrl_set(3:4, 3) = mag;
        t_end = 40; fsm.ctrl_set(end,1) = t_end; % See roll damping and drift
        
    case 'rudder_kick'
        mag = 5.0 * pi/180; % 5 deg left
        fsm.ctrl_set(3:4, 4) = mag;
        t_end = 60; fsm.ctrl_set(end,1) = t_end; % See Dutch Roll damping
        
    case 'flap_deployment'
        % Deploy 15 degrees of flaps at t=2
        fsm.delta_flap_static = 15 * pi/180; 
        t_end = 120; fsm.ctrl_set(end,1) = t_end; % See long-term trim shift
end

fprintf('Simulating scenario: %s...\n', scenario);
[tvec, Y] = ode45(@(t,y) simsub3d(t, y, fsm), [0 t_end], x3d(1:12));

% --- PLOT RESULTS ---
figure(1); clf; set(gcf, 'Color', 'w');

subplot(3,2,1); plot(tvec, Y(:,1)*3.6); grid on; ylabel('IAS (km/h)'); title('Airspeed');
subplot(3,2,2); plot(tvec, Y(:,12)); grid on; ylabel('Alt (m)'); title('Altitude');
subplot(3,2,3); plot(tvec, Y(:,7)*180/pi); grid on; ylabel('Bank (deg)'); title('Roll Angle (\phi)');
subplot(3,2,4); plot(tvec, Y(:,8)*180/pi); grid on; ylabel('Pitch (deg)'); title('Pitch Angle (\theta)');
subplot(3,2,5); plot(tvec, Y(:,9)*180/pi); grid on; ylabel('Heading (deg)'); title('Yaw Angle (\psi)');

subplot(3,2,6); hold on; grid on;
t_eval = linspace(0, t_end, 100);
ctrls = interp1(fsm.ctrl_set(:,1), fsm.ctrl_set(:,2:5), t_eval);
plot(t_eval, ctrls(:,1)*180/pi, 'r', 'LineWidth', 1.2);
plot(t_eval, ctrls(:,2)*180/pi, 'b', 'LineWidth', 1.2);
plot(t_eval, ctrls(:,3)*180/pi, 'g', 'LineWidth', 1.2);
legend('Elev', 'Ail', 'Rud', 'Location', 'best'); ylabel('Input (deg)'); title('Control History');
xlabel('Time (s)');
