% verify_downwash.m
% Extracts d_epsilon / d_alpha at the tail location in NLLT

addpath('../matlab');
acg = make_windex();

uoo = 30;
alpha_rad = [0, 2] * pi/180;
fs0 = flight_state(acg, uoo, 0, alpha_rad(1), 0);
fs1 = flight_state(acg, uoo, 0, alpha_rad(2), 0);

% Solve for wing lift distribution
fprintf('Solving for alpha=0...\n');
fs0.gamma = pointsolve(acg, fs0);
fprintf('Solving for alpha=2...\n');
fs1.gamma = pointsolve(acg, fs1);

% Tail location (Horizontal Tail AC approx)
x_tail = 2.8; 
y_tail = 0.0;
z_tail = 0.98;

% Get induced velocity at tail location
[~, v_t0] = segment_coefficients(acg, fs0, fs0.gamma);
[~, v_t1] = segment_coefficients(acg, fs1, fs1.gamma);

% Actually, segment_coefficients returns velocity for all segments.
% Let's find the tail segments (rt)
rt = acg.vxhtp;
vt0 = v_t0(rt, :);
vt1 = v_t1(rt, :);

% Average vertical velocity w (z-axis in vortex coords depends on convention)
% In segment_coefficients, vt is the total velocity vector in body axes?
% Body axes: x-back, y-right, z-up.
% Downwash epsilon = arctan(w / u).
w0 = mean(vt0(:, 3)); % z-component (upward)
u0 = mean(vt0(:, 1)); % x-component (forward)
eps0 = -atan2(w0, u0); % positive downwash

w1 = mean(vt1(:, 3));
u1 = mean(vt1(:, 1));
eps1 = -atan2(w1, u1);

deps_da = (eps1 - eps0) / (alpha_rad(2) - alpha_rad(1));

fprintf('\nDOWNWASH AT HT (z=0.98m):\n');
fprintf('  alpha=0: eps = %.4f deg\n', eps0 * 180/pi);
fprintf('  alpha=2: eps = %.4f deg\n', eps1 * 180/pi);
fprintf('  d_epsilon / d_alpha: %.4f\n', deps_da);

% Compare to wing plane (z=0)
z_tail_low = 0.1; % Near wing plane
% Need a point solver or velocity probe? 
% Let's use the induced velocity function directly if available.
% Or just look at the wing segments downwash?
rw = acg.vxwing;
% Induced at wing itself (mid-panel)
at0 = vortex_alpha(acg, v_t0);
at1 = vortex_alpha(acg, v_t1);
% This alpha includes farfield alpha.
% Ind alpha_ind = alpha - localpha?
% alpha_ind = d_epsilon
deps_da_wing = (at1(rw) - at0(rw)) / (alpha_rad(2) - alpha_rad(1)); % This is crude

% Exit
exit;
