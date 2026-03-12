% verify_windex_geometry.m
% Double checks the Windex geometry from the make_windex.m definition.
% Calculates Area, Span, and Mean Aerodynamic Chord (MAC).

clear all;
close all;
clc;

% Add paths
addpath('c:\Users\teodo\Desktop\flight_mechanics\nonlinear_lifting_line\matlab');

fprintf('=== WINDEX GEOMETRY VERIFICATION ===\n\n');

% Load geometry
acg = make_windex();

% Reference values from make_windex
S_ref_val = acg.Sref;
b_ref_val = acg.bref;
c_ref_val = acg.cref;

% Main wing indices
rw = acg.vxwing;

% Extract segment data for the wing
pa = acg.pa(rw,:);
pb = acg.pb(rw,:);
c_avg = acg.chord(rw);

% 1. Calculate Span (b)
y_all = [pa(:,2); pb(:,2)];
span = max(y_all) - min(y_all);

% 2. Calculate Area (S) and MAC
% We need to reconstruct the chords at the endpoints for better precision,
% but since add_segments stores c_avg = (c1+c2)/2, and area = c_avg * dy,
% the area sum is exact for trapezoids.

dy = sqrt(sum((pb - pa).^2, 2));
S_segs = c_avg .* dy;
S_total = sum(S_segs);

% For MAC (mean aerodynamic chord): c_bar = (1/S) * integral(c^2 dy)
% For a trapezoid: integral(c^2 dy) = dy/3 * (c1^2 + c1c2 + c2^2)
% We don't have c1 and c2 directly here, but we can approximate it well
% using c_avg if segments are small, or we can try to find the stations.
% Let's use c_avg as an approximation first: c_bar approx = sum(c_avg * S_seg) / S
c_bar_approx = sum(c_avg .* S_segs) / S_total;

% Accurate MAC LE calculation
% X_le_mac = (1/S) * integral(x_le(y) * c(y) dy)
% x_le at segment center is approx: (pa_x + pb_x)/2 - 0.25*c_avg 
% (reversing the 1/4 chord shift in add_segments)
x_le_seg = 0.5 * (pa(:,1) + pb(:,1)) - 0.25 * c_avg;
X_le_mac = sum(x_le_seg .* S_segs) / S_total;

% 3. Results Comparison
AR = span^2 / S_total;

fprintf('--- Calculated Values ---\n');
fprintf('Span (b):      %8.4f m\n', span);
fprintf('Area (S):      %8.4f m^2\n', S_total);
fprintf('AR:            %8.4f\n', AR);
fprintf('MAC (approx):  %8.4f m\n', c_bar_approx);
fprintf('MAC LE (X):    %8.4f m from Datum\n', X_le_mac);

fprintf('\n--- Reference Values in acg ---\n');
fprintf('Sref:          %8.4f m^2\n', S_ref_val);
fprintf('bref:          %8.4f m\n', b_ref_val);
fprintf('cref (MAC):    %8.4f m\n', c_ref_val);

fprintf('\n--- Discrepancies ---\n');
fprintf('Area error:    %8.4f %%\n', (S_total - S_ref_val)/S_ref_val * 100);
fprintf('Span error:    %8.4f %%\n', (span - b_ref_val)/b_ref_val * 100);
fprintf('MAC error:     %8.4f %%\n', (c_bar_approx - c_ref_val)/c_ref_val * 100);

if abs(S_total - S_ref_val)/S_ref_val > 0.05
    fprintf('\nWARNING: Significant area discrepancy detected!\n');
end
