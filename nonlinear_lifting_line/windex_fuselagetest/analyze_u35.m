% analyze_u35.m
% Simple script to extract CD vs alpha and Sref from u35.log

% Setup paths
p1 = '../../Windtunnel_fuse_w_tail/';
p2 = '../../Windtunnel_fuselagetest/';
addpath(p1);
log_file = fullfile(p2, 'u35.log');

fprintf('Searching for readwingloads.m in: %s\n', fullfile(pwd, p1));
if ~exist('readwingloads.m', 'file')
    fprintf('ERROR: readwingloads.m not found in path!\n');
end

fprintf('Searching for log file: %s\n', fullfile(pwd, log_file));
if ~exist(log_file, 'file')
    fprintf('ERROR: Log file not found!\n');
end

fprintf('Reading %s...\n', log_file);
[abq, ~, ~, coeff, span, Sref, cref] = readwingloads(log_file);

alpha = abq(:,1); % [deg]
CD = coeff(:,1);

fprintf('\n--- Wind Tunnel Log Info ---\n');
fprintf('Reference Area (Sref): %.4f m^2\n', Sref);
fprintf('Reference Chord (cref): %.4f m\n', cref);
fprintf('Reference Span:         %.4f m\n', span);

fprintf('\n--- CD vs Alpha Data ---\n');
fprintf('%10s | %10s\n', 'Alpha [deg]', 'CD');
fprintf('----------------------------\n');
for i = 1:length(alpha)
    fprintf('%10.2f | %10.5f\n', alpha(i), CD(i));
end

% Simple plot
figure('Color', 'w');
plot(alpha, CD, 'b-o', 'LineWidth', 1.5);
grid on; xlabel('\alpha [deg]'); ylabel('C_D');
title(['Experimental Drag Coefficient (u=35 m/s, Sref = ' num2str(Sref) ' m^2)']);
saveas(gcf, 'u35_drag_alpha.png');
fprintf('\nPlot saved as u35_drag_alpha.png\n');
