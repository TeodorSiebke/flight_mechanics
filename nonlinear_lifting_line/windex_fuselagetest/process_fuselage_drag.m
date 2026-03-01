function process_fuselage_drag()
    % Script to isolate fuselage drag from wind tunnel log files
    % by subtracting the NLL-calculated drag of the chopped wing.
    
    close all; clc;
    
    % --- Paths ---
    % Add the folder with readwingloads.m
    addpath('../../Windtunnel_fuse_w_tail/');
    % Add NLL library
    addpath('../matlab/');
    
    log_dir = '../../Windtunnel_fuselagetest/';
    log_files = {'u15.log', 'u25.log', 'u35.log', 'u40.log'};
    
    % --- Load Chopped Wing Model ---
    fprintf('Initializing Chopped Wing NLL Model...\n');
    acg = make_windexfusewing();
    
    % Results storage
    results = struct();
    
    for f = 1:length(log_files)
        fname = [log_dir log_files{f}];
        
        if ~exist(fname, 'file')
            %{
            % Alternative location check if pathing is relative to execution dir
            fname = ['../../windtunnel_fuselagetest/' log_files{f}];
            if ~exist(fname, 'file')
                continue;
            end
            %}
            continue;
        end
        
        fprintf('\nProcessing %s...\n', log_files{f});
        
        % Read log
        [abq, ~, ~, coeff, span, Sref_mod, cref_mod] = readwingloads(fname);
        
        alfa_exp = abq(:,1); % [deg]
        q_exp = abq(:,3);    % [Pa]
        CD_exp = coeff(:,1);
        CL_exp = coeff(:,3);
        
        % Determine approximate speed
        V_avg = mean(sqrt(q_exp / (0.5 * 1.225)));
        fprintf('  Average speed: %.1f m/s\n', V_avg);
        
        % --- Run NLL for chopped wing ---
        % Create flight state
        fs = flight_state(acg, V_avg, 0, 0, 0);
        fs.delta_elevator = 0; fs.delta_rudder = 0; fs.delta_flap = 0; fs.delta_aileron = 0;
        
        % Limit alfa to avoid convergence issues at high angles
        idx_valid = find(alfa_exp > -2 & alfa_exp < 2);
        alfa_sim = alfa_exp(idx_valid);
        
        % Run plainpolar
        [CL_sim_valid, CD_sim_valid, ~] = plainpolar(acg, fs, alfa_sim * pi/180);
        
        % Reconstruct full CD_sim with NaNs
        CD_sim = NaN(size(alfa_exp));
        CD_sim(idx_valid) = CD_sim_valid;
        
        % --- Fuselage Drag Calculation ---
        CD_fuse = CD_exp - CD_sim;
        
        % Find value at alfa ~ 0 (or average near cruise)
        idx_zero = find(abs(alfa_exp) < 0.5);
        if ~isempty(idx_zero)
            CD_fuse_zero = mean(CD_fuse(idx_zero));
            fprintf('  Fuselage CD (model scale) at alpha~0: %.5f\n', CD_fuse_zero);
            fprintf("CD_exp = %.5f\n", mean(CD_exp(idx_zero)));
            fprintf("CD_sim = %.5f\n", mean(CD_sim(idx_zero)));
            
            % Scale to full scale
            % S*CD_full = (CD_fuse_zero * Sref_mod) * 25
            SCD_full = (CD_fuse_zero * Sref_mod) * 25;
            fprintf('  Full scale S*CD contribution: %.5f m^2\n', SCD_full);
        end
        
        % Store for plotting
        results(f).V = V_avg;
        results(f).alfa = alfa_exp;
        results(f).CD_exp = CD_exp;
        results(f).CD_sim = CD_sim';
        results(f).CD_fuse = CD_fuse;
    end
    
    % --- Plotting ---
    if isempty(fields(results))
        fprintf('No data files processed. Check paths.\n');
        return;
    end
    
    figure('Color', 'w', 'Position', [100 100 800 600]);
    subplot(2,1,1);
    hold on; grid on;
    for f = 1:length(log_files)
        if f <= length(results) && ~isempty(results(f).alfa)
            plot(results(f).alfa, results(f).CD_exp, 'o-', 'DisplayName', ['Total ' log_files{f}]);
            plot(results(f).alfa, results(f).CD_sim, '--k', 'HandleVisibility', 'off');
        end
    end
    ylabel('C_D Total and Wing-Only');
    legend('Location', 'best');
    title('Wind Tunnel Drag (Total) vs Wing-Only NLL (Dashed)');
    
    subplot(2,1,2);
    hold on; grid on;
    for f = 1:length(log_files)
        if f <= length(results) && ~isempty(results(f).alfa)
            plot(results(f).alfa, results(f).CD_fuse, 'd-', 'DisplayName', ['Fuselage ' log_files{f}]);
        end
    end
    xlabel('\alpha [deg]');
    ylabel('C_D Fuselage');
    legend('Location', 'best');
    title('Extracted Fuselage Drag Coefficient');
    
    saveas(gcf, 'Fuselage_Drag_Analysis.png');
    fprintf('\nAnalysis complete. Plot saved.\n');
end
