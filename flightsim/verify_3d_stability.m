% verify_3d_stability.m
% Verification of 3D stability for landing approach: 
% -3 deg flight path, 25 m/s, 20 deg flaps.

clear; close all;
addpath(pwd);
fprintf('--- 3D Landing Stability Analysis (-3 deg Glide Path) ---\n');
disp('Starting stability analysis...');

% 1. Load model with 20 degree flaps
fsm = make_fsim(0);

% 2. Target Flight Condition
vset = 25; % m/s
gamma_target = -3 * pi/180; % Target flight path angle
fprintf('Airspeed: %.1f m/s\n', vset);
fprintf('Target Glide Slope: %.1f deg\n', gamma_target * 180/pi);

% CG setting (Nominal 9% SM)
x_np = 0.43;
cref = fsm.cref;
fsm.cog(1) = x_np - 0.09 * cref;

% Initial guess: [u v w p q r phi theta psi xed yed alt deltae deltaa deltar deltap]
% Note: deltaa, deltar usually 0 for symmetric trim
x = [vset 0.0 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1000 -0.1 0.0 0.0 0.2]';

% 3. Trim for Glide Slope
% To get -3 deg flight path, we need to satisfy: u_dot=0, w_dot=0, q_dot=0, but also 
% altdot = -3 deg path at vset
% altdot = vspeed * sin(-gamma) ? No, fplmod3d altdot is defined such that positive is UP?
% In fplmod3d: altdot = -tvbmx(3,1)*u - tvbmx(3,2)*v - tvbmx(3,3)*w
% For wings level, altdot = u*sin(theta) - w*cos(theta)
% Flight path angle gamma = atan2(-altdot, gspd) where gspd is horizontal.
% Or more simply: gamma = theta - alpha (for small angles).

% Trim variables: alpha, theta, delta_e, thrust (IVAR)
% Constraints: udot, wdot, qdot, gamma - target_gamma (IFUN)
ivar = [3 8 13 16]; % w, theta, delta_e, deltap
ifun = [1 3 5 12];  % udot, wdot, qdot, altdot (which relates to gamma)

x(1) = vset;
xtrim = x(ivar);

for iter = 1:100
    x(ivar) = xtrim;
    [xdot] = fplmod3d(0, x, fsm);
    
    % Target altdot for -3 deg path: altdot = vset * sin(gamma_target)
    target_altdot = vset * sin(gamma_target);
    
    ftrim = [xdot(1); xdot(3); xdot(5); xdot(12) - target_altdot];
    
    if norm(ftrim) < 1e-10, break; end
    
    % Jacobian
    J = zeros(4, 4);
    xstep = 1e-7;
    for j = 1:4
        xh = x; xh(ivar(j)) = xh(ivar(j)) + xstep;
        xmh = x; xmh(ivar(j)) = xmh(ivar(j)) - xstep;
        [xdoth] = fplmod3d(0, xh, fsm);
        [xdotmh] = fplmod3d(0, xmh, fsm);
        J(:,j) = ([xdoth(1); xdoth(3); xdoth(5); xdoth(12)] - ...
                  [xdotmh(1); xdotmh(3); xdotmh(5); xdotmh(12)]) / (2*xstep);
    end
    
    xtrim = xtrim - 0.5 * (J \ ftrim);
end
x(ivar) = xtrim;

% Check Results
alpha = atan2(x(3), x(1));
gamma_actual = x(8) - alpha;
fprintf('\nTrim Results:\n');
fprintf('  Alpha:      %.2f deg\n', alpha * 180/pi);
fprintf('  Theta:      %.2f deg\n', x(8) * 180/pi);
fprintf('  Gamma:      %.2f deg\n', gamma_actual * 180/pi);
fprintf('  Elevator:   %.2f deg\n', x(13) * 180/pi);
fprintf('  Thrust:     %.1f N (%.1f%%)\n', x(16)*fsm.thrustmax, x(16)*100);

if abs(gamma_actual - gamma_target) > 1e-4
    fprintf('WARNING: Glide slope trim failed to reach target precisely.\n');
end

% 4. Full 3D Stability Analysis
% State: [u v w p q r phi theta psi] (9 states)
st_idx = 1:9;
A = zeros(9, 9);
xstep = 1e-5;
for j = 1:9
    idx = st_idx(j);
    xh = x; xh(idx) = xh(idx) + xstep;
    xmh = x; xmh(idx) = xmh(idx) - xstep;
    [xdoth] = fplmod3d(0, xh, fsm);
    [xdotmh] = fplmod3d(0, xmh, fsm);
    A(:,j) = (xdoth(st_idx) - xdotmh(st_idx)) / (2*xstep);
end

eg = eig(A);

% 5. Separate Longitudinal and Lateral Modes
% Longitudinal: u, w, q, theta (states 1, 3, 5, 8)
Alon = A([1 3 5 8], [1 3 5 8]);
eg_lon = eig(Alon);

% Lateral: v, p, r, phi, psi (states 2, 4, 6, 7, 9)
Alat = A([2 4 6 7 9], [2 4 6 7 9]);
eg_lat = eig(Alat);

% Sort and display
fprintf('\n--- Longitudinal Modes ---\n');
[~, s_idx] = sort(abs(real(eg_lon)), 'descend');
eg_lon = eg_lon(s_idx);
for i=1:length(eg_lon)
    re = real(eg_lon(i));
    im = imag(eg_lon(i));
    fprintf('  Mode %d: %6.4f + %6.4fi', i, re, im);
    if re ~= 0
        t_val = log(2) / abs(re);
        if re < 0
            fprintf('  (T1/2 = %.2f s)', t_val);
        else
            fprintf('  (T2   = %.2f s)', t_val);
        end
    end
    fprintf('\n');
end

fprintf('\n--- Lateral Modes ---\n');
[~, s_idx] = sort(abs(real(eg_lat)), 'descend');
eg_lat = eg_lat(s_idx);
for i=1:length(eg_lat)
    re = real(eg_lat(i));
    im = imag(eg_lat(i));
    fprintf('  Mode %d: %6.4f + %6.4fi', i, re, im);
    if re ~= 0
        t_val = log(2) / abs(re);
        if re < 0
            fprintf('  (T1/2 = %.2f s)', t_val);
        else
            fprintf('  (T2   = %.2f s)', t_val);
        end
    end
    fprintf('\n');
end

% Mode Identification (Heuristic)
fprintf('\nMode Identification Summary:\n');
% Phugoid: lowest frequency complex in lon
[~, ph_idx] = min(abs(eg_lon)); % Simplification
fprintf('  Phugoid Eig:   %.4f + %.4fi\n', real(eg_lon(ph_idx)), imag(eg_lon(ph_idx)));

% Dutch Roll: complex pair in lat
dr_idx = find(imag(eg_lat) ~= 0);
if ~isempty(dr_idx)
    fprintf('  Dutch Roll Eig: %.4f + %.4fi\n', real(eg_lat(dr_idx(1))), imag(eg_lat(dr_idx(1))));
end

% Roll Subsidence: most stable real in lat
roll_idx = find(imag(eg_lat) == 0);
[~, rs_sub] = min(real(eg_lat(roll_idx))); 
fprintf('  Roll Sub Eig:  %.4f\n', eg_lat(roll_idx(rs_sub)));

% Spiral: least stable real in lat
[~, sp_sub] = max(real(eg_lat(roll_idx)));
fprintf('  Spiral Eig:    %.4f\n', eg_lat(roll_idx(sp_sub)));

% Stability Summary
if all(real(eg) < 0)
    fprintf('\n3D AIRCRAFT STATUS: STABLE\n');
else
    fprintf('\n3D AIRCRAFT STATUS: UNSTABLE (Check Spiral/Phugoid)\n');
end

% Plot Root Locus
figure('Color', 'w', 'Position', [100 100 800 600]);
plot(real(eg_lon), imag(eg_lon), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Longitudinal');
hold on; grid on;
plot(real(eg_lat), imag(eg_lat), 'bx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Lateral');
xlabel('Real Axis (1/s)'); ylabel('Imaginary Axis (rad/s)');
title('3D Stability Modes (Landing @ 25 m/s, -3 deg Path)');
legend('Location', 'best');
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
saveas(gcf, 'full_3d_stability_approach.png');

fprintf('\nAnalysis complete. Plot saved to full_3d_stability_approach.png\n');
