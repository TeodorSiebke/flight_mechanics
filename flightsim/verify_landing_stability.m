% verify_landing_stability.m
% Verification of the 20° flap aerodynamic database and stability at 25 m/s

clear; close all;
fprintf('--- Verifying Landing Stability (20 deg Flaps) ---\n');

% Load model with 20 degree flaps
fsm = make_fsim(20);

% Flight condition
vset = 25; % m/s
fprintf('Airspeed: %.1f m/s\n', vset);

% CG setting (Nominal 9% SM as in task6_analysis.m)
x_np = 0.43;
cref = fsm.cref;
fsm.cog(1) = x_np - 0.09 * cref;

% Initial guess: [u w q theta dist alt fuel de dp]
x = [vset 0.5 0.0 0.1 0 1000 0 -0.1 0.2]'; 

% Trim Loop
ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8];
xtrim = x(ivar);
for iter = 1:100 
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
    if norm(ftrim) < 1e-10, break; end
    xtrim = xtrim - 0.5 * (J \ ftrim); % damped newton
end
x(ivar) = xtrim;

fprintf('\nTrim State:\n');
fprintf('  Alpha:    %.2f deg\n', x(2)/x(1) * 180/pi);
fprintf('  Elevator: %.2f deg\n', x(8) * 180/pi);
fprintf('  Thrust:   %.1f N\n', x(9) * fsm.thrustmax);
fprintf('  Pitch:    %.2f deg\n', x(4) * 180/pi);

if norm(ftrim) > 1e-3
    fprintf('WARNING: Trim failed, residual norm: %.2e\n', norm(ftrim));
end

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
[~, s_idx] = sort(abs(real(eg)), 'descend');
eg = eg(s_idx);

fprintf('\nEigenvalues:\n');
for i = 1:4
    fprintf('  %.4f + %.4fi\n', real(eg(i)), imag(eg(i)));
end

% mode info
t_vals = log(2) ./ real(eg);
fprintf('\nMode Dynamics:\n');
fprintf('  Short Period T1/2: %.2f s\n', mean(t_vals(1:2)));
fprintf('  Phugoid T1/2 (or T2): %.2f s\n', mean(t_vals(3:4)));

% Stability check
if all(real(eg) < 0)
    fprintf('\nSTATUS: Stable\n');
else
    fprintf('\nSTATUS: Unstable\n');
end

% Plotting
figure('Color', 'w');
plot(real(eg), imag(eg), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
grid on; hold on;
xlabel('Real Axis (1/s)'); ylabel('Imaginary Axis (rad/s)');
title(sprintf('Landing Stability (20 deg Flaps) at %.0f m/s', vset));
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
saveas(gcf, 'landing_stability_verification.png');

fprintf('\nVerification complete. Plot saved to landing_stability_verification.png\n');
