% task6_step_response.m
% Longitudinal Step Response Simulation for the Windex
clear; close all;
addpath('flightsim');
fsm = make_fsim();

% --- CONFIGURATION: 9% Static Margin ---
x_np = 0.43;
cref = fsm.cref;
sm_target = 0.09;
xcg = x_np - sm_target * cref;
fsm.cog(1) = xcg;
fprintf('Simulating with 9%% SM (x_cg = %.3f m)\n', xcg);

% 1. Trim at 30 m/s level flight
vset = 30;
x_lon = [vset 0.4 0.0 0.05 0 1000 0 0.0 0.1]'; 
ivar_lon = [1 2 4 8 9]; ifun_lon = [1 2 3 6 8];
xtrim = x_lon(ivar_lon);
for iter = 1:50 
    x_lon(ivar_lon) = xtrim;
    [xdot_lon] = fplmod(0, x_lon, fsm);
    xstep = 1e-7; J = zeros(5, 5);
    for j = 1:5
        xh = x_lon; xh(ivar_lon(j)) = xh(ivar_lon(j)) + xstep;
        xmh = x_lon; xmh(ivar_lon(j)) = xmh(ivar_lon(j)) - xstep;
        [xdoth] = fplmod(0, xh, fsm); [xdotmh] = fplmod(0, xmh, fsm);
        J(:,j) = (xdoth(ifun_lon) - xdotmh(ifun_lon)) / (2*xstep);
    end
    ftrim = xdot_lon(ifun_lon); ftrim(5) = ftrim(5) - vset;
    if norm(ftrim) < 1e-9, break; end
    xtrim = xtrim - 0.7 * (J \ ftrim);
end
x_lon(ivar_lon) = xtrim;

% 2. Derive Linearized A and B matrices
% x = [u, w, q, theta, alt]
nx = 5; nu = 2;
idx_x = [1 2 3 4 6]; % u, w, q, theta, alt
idx_u = [8 9];       % de, dp
A = zeros(nx, nx); B = zeros(nx, nu);
xstep = 1e-5;

for j = 1:nx
    xh = x_lon; xh(idx_x(j)) = xh(idx_x(j)) + xstep;
    xmh = x_lon; xmh(idx_x(j)) = xmh(idx_x(j)) - xstep;
    [xdoth] = fplmod(0, xh, fsm); [xdotmh] = fplmod(0, xmh, fsm);
    A(:,j) = (xdoth(idx_x) - xdotmh(idx_x)) / (2*xstep);
end

for j = 1:nu
    xh = x_lon; xh(idx_u(j)) = xh(idx_u(j)) + xstep;
    xmh = x_lon; xmh(idx_u(j)) = xmh(idx_u(j)) - xstep;
    [xdoth] = fplmod(0, xh, fsm); [xdotmh] = fplmod(0, xmh, fsm);
    B(:,j) = (xdoth(idx_x) - xdotmh(idx_x)) / (2*xstep);
end

% 3. Simulation using ode45 (No toolbox required)
t_end = 600; % Extended to see altitude trend
t_span = [0 t_end];
h0 = x_lon(6);

% Scenario 1: Elevator Step (+1 deg = nose down nudge)
de_step = pi/180;
du_elev = [de_step; 0];
ode_elev = @(t, x) A*x + B*du_elev;
[t_elev, y_elev] = ode45(ode_elev, t_span, zeros(nx,1));
u_elev = vset + y_elev(:,1);
% Asymptote for Speed (using 4x4 part as 5th is integrator)
A4 = A(1:4, 1:4); B4 = B(1:4, :);
ss_elev = -A4 \ (B4(:,1) * de_step); 
u_ss_elev = vset + ss_elev(1);

% Scenario 2: Throttle Step (+0.1)
dp_step = 0.1;
du_throt = [0; dp_step];
ode_throt = @(t, x) A*x + B*du_throt;
[t_throt, y_throt] = ode45(ode_throt, linspace(0, t_end, 1000), zeros(nx,1));
h_throt = h0 + y_throt(:,5);
% For altitude, the asymptote is a ramp. We show the initial level and trend.
% Actually, let's see if the 5x5 is stable enough to give a height asymptote.
% (Due to density variation included in fplmod)
if cond(A) < 1e12
    ss_full = -A \ (B(:,2) * dp_step);
    h_ss_throt = h0 + ss_full(5);
else
    h_ss_throt = NaN;
end

% --- PLOTTING ---
figure(5); clf; set(gcf, 'Position', [100 100 900 700], 'Color', 'w');

subplot(2,1,1); hold on; grid on;
plot(t_elev, u_elev, 'b', 'LineWidth', 1.5, 'DisplayName', 'Response (Elevator +1^\circ)');
line([0 t_end], [u_ss_elev u_ss_elev], 'Color', 'r', 'LineStyle', '--', 'DisplayName', 'Speed Asymptote');
title('Airspeed Response to Elevator Step (+1^\circ)');
ylabel('Airspeed (m/s)'); legend('Location', 'best');

subplot(2,1,2); hold on; grid on;
plot(t_throt, h_throt, 'k', 'LineWidth', 1.5, 'DisplayName', 'Response (Throttle +0.1)');
if ~isnan(h_ss_throt)
    line([0 t_end], [h_ss_throt h_ss_throt], 'Color', 'r', 'LineStyle', '--', 'DisplayName', 'Height Asymptote');
end
title('Altitude Response to Throttle Step (+0.1)');
xlabel('Time (s)'); ylabel('Altitude (m)'); legend('Location', 'best');

sgtitle('Longitudinal Velocity and Altitude Step Response');
saveas(5, 'task6_step_response.png');
fprintf('Step response analysis complete.\n');
