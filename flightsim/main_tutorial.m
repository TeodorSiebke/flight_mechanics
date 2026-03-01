% main_tutorial.m
% A copy of the original flightsim template script, tweaked for Windex gliding.

fsm = make_fsim();

% --- 1. SETTINGS FOR GLIDING TRIM ---
% Trim variables: u, w, theta, elevator
ivar = [1 2 4 8]; 
% Trim functions: udot, wdot, qdot, airspeed
ifun = [1 2 3 8]; 

vset = 30.0; % Target airspeed (m/s)

% --- 2. INITIAL GUESS ---
% x: [u   w   q   theta  dist  alt   fuel  de  dp]
x = [vset 1.0 0.0 -0.05 0 1000 0 0 0]'; 

% --- 3. TRIM LOOP (Newton-Raphson) ---
xtrim = x(ivar);
fprintf('Trimming for V = %.1f m/s...\n', vset);

for iter = 1:20
    x(ivar) = xtrim;
    [xdot] = fplmod(0, x, fsm);
    
    xstep = 1e-6;
    J = zeros(length(ifun), length(ivar));
    for j = 1:length(ivar)
        xh = x; xh(ivar(j)) = xh(ivar(j)) + xstep;
        xmh = x; xmh(ivar(j)) = xmh(ivar(j)) - xstep;
        [xdoth] = fplmod(0, xh, fsm);
        [xdotmh] = fplmod(0, xmh, fsm);
        J(:,j) = (xdoth(ifun) - xdotmh(ifun)) / (2*xstep);
    end
    
    ftrim = xdot(ifun);
    ftrim(4) = ftrim(4) - vset; 
    
    % Use a damped step for better stability
    xtrim = xtrim - 0.5 * (J \ ftrim);
    
    fprintf('  Iteration %d: error = %e\n', iter, norm(ftrim));
    if norm(ftrim) < 1e-8, break; end
end

x(ivar) = xtrim;
[xdot] = fplmod(0, x, fsm); % Final call to get steady state
fprintf('Trim Complete: Alpha = %.2f deg, DE = %.2f deg\n', xdot(9)*180/pi, x(8)*180/pi);

% --- 4. FLIGHT SIMULATION (Time Domain) ---
fsm.dp_set = x(9); % Set throttle setting for simulation
% Define a control input history (Step input to elevator)
% Column 1: Time (s), Column 2: Elevator deflection (rad)
fsm.de_set = [ 0   x(8)
               10  x(8)
               11  x(8) - 2*pi/180  % Nose-up nudge at t=10s
               20  x(8) - 0.5*pi/180
               21  x(8)               % Return to trim at t=20s
             1000  x(8)];

fprintf('Running time simulation for 60 seconds...\n');
[tvec, Y] = ode45(@(t,y) simsub(t, y, fsm), [0 60], x(1:7));

% --- 5. PLOTTING ---
figure(1);
subplot(3,1,1); plot(tvec, Y(:,1)); grid on; ylabel('u (m/s)'); title('Windex Response to Elevator Nudge');
subplot(3,1,2); plot(tvec, Y(:,2)); grid on; ylabel('w (m/s)'); 
subplot(3,1,3); plot(tvec, Y(:,4)*180/pi); grid on; ylabel('theta (deg)'); xlabel('Time (s)');

figure(2);
plot(Y(:,5), Y(:,6)); grid on; xlabel('Distance (m)'); ylabel('Altitude (m)');
title('Gliding Trajectory'); axis equal;
