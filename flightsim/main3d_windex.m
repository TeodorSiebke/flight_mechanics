% main3d_windex.m
% Full 3D simulation of the Windex aircraft using 12-state rigid body equations.

fsm = make_fsim();

% --- 1. SETTINGS FOR 3D TRIM ---
% Trim at V = 30 m/s, Level flight.
vset = 30.0;

% x3d: [u v w p q r phi theta psi xe ye alt de da dr dp]
% Indices: 1 2 3 4 5 6   7     8   9 10 11  12 13 14 15 16
x3d = zeros(16,1);
x3d(1) = vset;
x3d(12) = 1000; % Altitude 1000m

% Select variables and functions for trim
% ivar: [1: u, 3: w, 8: theta, 13: de, 16: dp] (Longitudinal trim)
ivar = [1 3 8 13 16]; 
ifun = [1 3 5 8 12]; % [udot, wdot, qdot, alt_dot (climb), airspeed]

xtrim = [vset 1.0 -0.05 0.0 0.1]'; % Initial guess

fprintf('=== Windex 3D Trim Analysis ===\n');
for iter = 1:15
    x3d(ivar) = xtrim;
    [xdot] = fplmod3d(0, x3d, fsm);
    vcurr = sqrt(x3d(1)^2 + x3d(2)^2 + x3d(3)^2);
    
    xstep = 1e-6;
    J = zeros(5, 5);
    for j = 1:5
        xh = x3d; xh(ivar(j)) = xh(ivar(j)) + xstep;
        xmh = x3d; xmh(ivar(j)) = xmh(ivar(j)) - xstep;
        
        [xdoth] = fplmod3d(0, xh, fsm);
        [xdotmh] = fplmod3d(0, xmh, fsm);
        
        vhr = sqrt(xh(1)^2 + xh(2)^2 + xh(3)^2);
        vmr = sqrt(xmh(1)^2 + xmh(2)^2 + xmh(3)^2);
        
        % Combine state derivatives and airspeed derivative
        der_h = [xdoth(1); xdoth(3); xdoth(5); xdoth(12); vhr];
        der_m = [xdotmh(1); xdotmh(3); xdotmh(5); xdotmh(12); vmr];
        J(:,j) = (der_h - der_m) / (2*xstep);
    end
    
    % Residual vector: [udot, wdot, qdot, altdot, airspeed-vset]
    ftrim = [xdot(1); xdot(3); xdot(5); xdot(12); vcurr-vset];
    
    xtrim = xtrim - 0.5 * (J \ ftrim); 
    if norm(ftrim) < 1e-8, break; end
end

x3d(ivar) = xtrim;
fprintf('Trim Results (3D Model):\n');
fprintf('  Speed: %.1f m/s, Alpha: %.2f deg, Theta: %.2f deg\n', vset, atan(x3d(3)/x3d(1))*180/pi, x3d(8)*180/pi);
fprintf('  Elevator: %.2f deg, Throttle: %.2f\n', x3d(13)*180/pi, x3d(16));

% --- 2. DYNAMIC MANEUVER: STEP TURN ---
% Define control history for 3D simulation
% tvec, [de, da, dr, dp]
fsm.ctrl_set = [ 0   x3d(13)  0.0       0.0  x3d(16)
                 5   x3d(13)  0.0       0.0  x3d(16)
                 6   x3d(13)  2*pi/180  0.0  x3d(16) % Roll command (Stick right)
                 10  x3d(13)  2*pi/180  0.0  x3d(16)
                 11  x3d(13)  0.0       0.0  x3d(16) % Return stick to neutral
                 200 x3d(13)  0.0       0.0  x3d(16) ];

fprintf('\nRunning 3D time simulation (Roll Command)...\n');
[tvec, Y] = ode45(@(t,y) simsub3d(t, y, fsm), [0 60], x3d(1:12));

% --- 3. PLOTTING ---
figure(1);
subplot(3,1,1); plot(tvec, Y(:,1)); grid on; ylabel('u (m/s)'); title('3D Response: Aileron Pulse');
subplot(3,1,2); plot(tvec, Y(:,7)*180/pi); grid on; ylabel('Bank Angel Phi (deg)');
subplot(3,1,3); plot(tvec, Y(:,9)*180/pi); grid on; ylabel('Heading Psi (deg)');
xlabel('Time (s)');

figure(2);
plot(Y(:,11), Y(:,10)); grid on; xlabel('East (m)'); ylabel('North (m)');
title('Ground Track during Turn'); axis equal;
