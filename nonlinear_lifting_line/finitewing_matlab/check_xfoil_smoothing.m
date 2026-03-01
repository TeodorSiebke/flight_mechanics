%
% Diagnostic script to check the smoothness of XFOIL polars
%

clear all;
close all;

% --- Configuration ---
% Add core NLLLT library for importing
addpath('../matlab');
polar_dir = '../windex_xfoil/';

% Check a wide range of Reynolds numbers
max_re = 20;
cmap = lines(max_re);

for i = 1:max_re
    fname = sprintf('%sh15_re%de5.txt', polar_dir, i);
    
    if exist(fname, 'file') == 2
        try
            [alpha, CL, CD, CM] = import_xfoil_polar(fname);
            alpha_deg = alpha * 180/pi;
            
            % Plot CL
            figure(1);
            plot(alpha_deg, CL, '.-', 'Color', cmap(i,:), 'DisplayName', sprintf('Re=%de5', i));
            hold on;
            
            % Plot CD
            figure(2);
            plot(alpha_deg, CD, '.-', 'Color', cmap(i,:), 'DisplayName', sprintf('Re=%de5', i));
            hold on;
            
            % Plot CM
            figure(3);
            plot(alpha_deg, CM, '.-', 'Color', cmap(i,:), 'DisplayName', sprintf('Re=%de5', i));
            hold on;
            
        catch ME
            fprintf('Error reading %s: %s\n', fname, ME.message);
        end
    else
        fprintf('File not found: %s\n', fname);
    end
end

% Finalize plots
figure(1);
xlabel('Angle of Attack \alpha [deg]');
ylabel('Lift Coefficient C_L');
title('XFOIL Smoothing Check: C_L vs \alpha');
grid on; legend('show', 'Location', 'SouthEast');

figure(2);
xlabel('Angle of Attack \alpha [deg]');
ylabel('Drag Coefficient C_D');
title('XFOIL Smoothing Check: C_D vs \alpha');
grid on; legend('show', 'Location', 'NorthWest');

figure(3);
xlabel('Angle of Attack \alpha [deg]');
ylabel('Moment Coefficient C_m');
title('XFOIL Smoothing Check: C_m vs \alpha');
grid on; legend('show', 'Location', 'SouthWest');

fprintf('Check plots generated. Look for sharp "kinks" or jumps in the curves.\n');
