%% Task 1: Geometry Analysis (MAC Calculation)
% Calculates Mean Aerodynamic Chord (MAC) and its leading edge position.
clear; clc;

% Load Configuration
if ~exist('Config', 'var')
    windex_config;
end

fprintf('\n=== TASK 1: GEOMETRY ANALYSIS ===\n');

% CALCULATIONS
[MAC_W, XLE_W] = calculate_MAC(Config.Geom.Wing);
[MAC_H, XLE_H] = calculate_MAC(Config.Geom.Tail);

% OUTPUT RESULTS
fprintf('----------------------------------------------------\n');
fprintf('SURFACE          MAC LENGTH (m)    X-LE OF MAC (m)  \n');
fprintf('----------------------------------------------------\n');
fprintf('Main Wing:         %8.4f          %8.4f\n', MAC_W, XLE_W);
fprintf('Horiz. Tail:       %8.4f          %8.4f\n', MAC_H, XLE_H);
fprintf('----------------------------------------------------\n');
fprintf('Note: X-LE is the distance aft of that surface''s root LE.\n');

% HELPER FUNCTION
function [MAC, XLE] = calculate_MAC(sections)
    total_area = 0;
    weighted_mac = 0;
    weighted_xle = 0;

    for i = 1:(size(sections, 1) - 1)
        y1 = sections(i, 1);    y2 = sections(i+1, 1);
        c1 = sections(i, 2);    c2 = sections(i+1, 2);
        x1 = sections(i, 3);    x2 = sections(i+1, 3);
        
        % Segment properties
        b_seg = y2 - y1;
        taper = c2 / c1;
        area_seg = b_seg * (c1 + c2) / 2;
        
        % Local MAC and its X-LE location for this trapezoid
        mac_seg = (2/3) * c1 * (1 + taper + taper^2) / (1 + taper);
        xle_seg = x1 + (x2 - x1) * (c1 + 2*c2) / (3 * (c1 + c2));
        
        % Area-weighting
        total_area = total_area + area_seg;
        weighted_mac = weighted_mac + (mac_seg * area_seg);
        weighted_xle = weighted_xle + (xle_seg * area_seg);
    end
    
    MAC = weighted_mac / total_area;
    XLE = weighted_xle / total_area;
end
