% Main gliding simulation and polar sweep
fsm = make_fsim();

% Select variables and function for GLIDING trim
% ivar: [1: u, 2: w, 4: theta, 8: de]
% ifun: [1: udot, 2: wdot, 3: qdot, 8: airspd]
ivar = [1 2 4 8];
ifun = [1 2 3 8];

v_sweep = 23:1:65; % Start above stall speed (~21-22 m/s)
sink_rates = zeros(length(v_sweep), 1);
alphas = zeros(length(v_sweep), 1);
elevators = zeros(length(v_sweep), 1);

% Initial guess for the first point (High speed is usually easier to trim)
% x: [u w q theta dist alt fuel de dp]
x = [v_sweep(1) 1.0 0.0 -0.05 0 1000 0 0 0]'; 

fprintf('=== Windex Gliding Polar Sweep (Flight Sim) ===\n');
fprintf('  V (m/s) |  Sink (m/s) | Alpha (deg) | DE (deg)\n');
fprintf('-------------------------------------------------\n');

for i = 1:length(v_sweep)
    vset = v_sweep(i);
    xtrim = x(ivar);
    x(9) = 0; % Throttle set to zero for glide
    
    % Newton iteration for trim
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
        ftrim(4) = ftrim(4) - vset; % Constraint: airspeed = vset
        
        % Damped Newton step to prevent "wild" divergence
        delta = J \ ftrim;
        xtrim = xtrim - 0.5 * delta; 
        
        if norm(ftrim) < 1e-7, break; end
    end
    
    % Store results
    [xdot_final] = fplmod(0, x, fsm);
    sink_rates(i) = -xdot_final(6);
    alphas(i) = xdot_final(9) * 180/pi;
    elevators(i) = x(8) * 180/pi;
    
    fprintf('%9.1f | %10.4f | %11.2f | %8.2f\n', vset, sink_rates(i), alphas(i), elevators(i));
end

% Plotting
figure(1);
plot(v_sweep, sink_rates, 'k-o', 'LineWidth', 2);
grid on; xlabel('Airspeed (m/s)'); ylabel('Sink Rate (m/s)');
title('Windex Gliding Polar (flight sim)');
set(gca, 'YDir', 'reverse');

figure(2);
subplot(2,1,1); plot(v_sweep, alphas, 'r', 'LineWidth', 2); grid on; ylabel('Alpha (deg)');
subplot(2,1,2); plot(v_sweep, elevators, 'g', 'LineWidth', 2); grid on; ylabel('Elevator (deg)');
xlabel('Airspeed (m/s)');
