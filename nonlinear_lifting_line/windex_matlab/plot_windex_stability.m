% plot_windex_stability.m
% Copy of geometry plot script with Aerodynamic Centers and Neutral Point added.
% Based on plot_windex.m and calculate_total_np_nll.m

clear all;
close all;

%% 1. GEOMETRY DATA (from plot_windex.m)
% Patches for flat panel model of the windex, twist not shown
% Each row gives x1 y1 z1 chord1 x2 y2 z2 chord2 for each patch 
pxyz=[
 0.1600   -6.0000    0.3140    0.2342    0.1300   -5.5000    0.2878    0.3000
 0.1300   -5.5000    0.2878    0.3000    0.0500   -3.5000    0.1832    0.5233
 0.0500   -3.5000    0.1832    0.4933         0   -0.3000    0.0157    0.5800
      0   -0.3000    0.0157    0.7283         0         0         0    0.7283
      0         0         0    0.7283         0    0.3000    0.0157    0.7283
      0    0.3000    0.0157    0.5800    0.0500    3.5000    0.1832    0.4933
 0.0500    3.5000    0.1832    0.5233    0.1300    5.5000    0.2878    0.3000
 0.1300    5.5000    0.2878    0.3000    0.1600    6.0000    0.3140    0.2342
 0.3942   -6.0000    0.3140    0.0658    0.4300   -5.5000    0.2878    0.1000
 0.4300   -5.5000    0.2878    0.1000    0.5733   -3.5000    0.1832    0.1182
 0.5433   -3.5000    0.1832    0.1483    0.5800   -0.3000    0.0157    0.1483
 0.5800    0.3000    0.0157    0.1483    0.5433    3.5000    0.1832    0.1483
 0.5733    3.5000    0.1832    0.1182    0.4300    5.5000    0.2878    0.1000
 0.4300    5.5000    0.2878    0.1000    0.3942    6.0000    0.3140    0.0658
 2.9215   -1.0750    0.9800    0.1865    2.7325         0    0.9800    0.3755
 2.7325         0    0.9800    0.3755    2.9215    1.0750    0.9800    0.1865
 3.1080   -1.0750    0.9800    0.0785    3.1080         0    0.9800    0.1595
 3.1080         0    0.9800    0.1595    3.1080    1.0750    0.9800    0.0785
 1.6314         0    0.0120    0.8413    2.7325         0    0.9800    0.3500
 2.4727         0    0.0120    0.3605    3.0825         0    0.9800    0.1500];

%% 2. STABILITY DATA (from calculate_total_np_nll.m results)
% Positions in meters from Wing Root LE (Datum x=0)
Xac_w = 0.2185; 
Xac_t = 2.9271;
X_np  = 0.427;

%% 3. PLOTTING
figure('Name', 'Windex Geometry and Stability Points', 'Color', 'w');
hold on;
grid on;
axis equal;
xlabel('x [m] (Aft)');
ylabel('y [m] (Right)');
zlabel('z [m] (Up)');
title('Windex Geometry: AC Positions and Neutral Point');

% Plot Aircraft Geometry
for ip=1:size(pxyz,1)
    xp=[pxyz(ip,1) pxyz(ip,5)  pxyz(ip,5)+pxyz(ip,8)  pxyz(ip,1)+pxyz(ip,4)];
    yp=[pxyz(ip,2) pxyz(ip,6)  pxyz(ip,6) pxyz(ip,2) ];
    zp=[pxyz(ip,3) pxyz(ip,7)  pxyz(ip,7) pxyz(ip,3) ];
    xp(5)=xp(1); yp(5)=yp(1); zp(5)=zp(1);
    plot3(xp, yp, zp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end

% Plot Stability Markers
% We place them on the symmetry plane (y=0) for visibility
h1 = plot3(Xac_w, 0, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Wing AC');
h2 = plot3(Xac_t, 0, 0.98, 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Tail AC');
h3 = plot3(X_np,  0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Total NP');

% Add vertical lines to show X-position clearly
line([Xac_w Xac_w], [0 0], [-0.5 0.5], 'Color', 'r', 'LineStyle', '--', 'HandleVisibility', 'off');
line([Xac_t Xac_t], [0 0], [0.5 1.5], 'Color', 'b', 'LineStyle', '--', 'HandleVisibility', 'off');
line([X_np  X_np],  [0 0], [-0.5 0.5], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

legend([h1 h2 h3], 'Location', 'best');
view(45, 30);

fprintf('Stability Points (from Root LE):\n');
fprintf('Wing AC: %.4f m\n', Xac_w);
fprintf('Tail AC: %.4f m\n', Xac_t);
fprintf('Neutral Point: %.4f m\n', X_np);
