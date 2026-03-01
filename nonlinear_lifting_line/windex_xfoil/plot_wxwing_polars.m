function plot_wxwing_polars()
    % PLOT_WXWING_POLARS Plots CL and CM vs Alpha for Windex wing XFOIL results
    % for all found deflections and Reynolds numbers.

    % Get the directory of the script
    script_dir = fileparts(mfilename('fullpath'));
    
    % Find all wxwing*.txt files
    files = dir(fullfile(script_dir, 'wxwing*.txt'));
    
    if isempty(files)
        fprintf('No wxwing*.txt files found in %s\n', script_dir);
        return;
    end

    % Parse files to extract Re and Deflection
    data_struct = struct();
    for i = 1:length(files)
        fname = files(i).name;
        
        % Extract Re (e.g., re15e5 or re3e5)
        re_match = regexp(fname, 're(\d+e\d+)', 'tokens');
        if isempty(re_match)
            continue;
        end
        re_str = re_match{1}{1};
        
        % Extract Deflection (e.g., d10 or d-20)
        d_match = regexp(fname, '_d(-?\d+)_', 'tokens');
        if isempty(d_match)
            deflection = 0; % Default if no _dXX_ in name
        else
            deflection = str2double(d_match{1}{1});
        end
        
        % Read data
        file_path = fullfile(script_dir, fname);
        try
            raw = importdata(file_path, ' ', 12);
            if isstruct(raw)
                d = raw.data;
            else
                d = readmatrix(file_path, 'NumHeaderLines', 12);
            end
            
            % Store in structured format: data_struct.(re_str).(deflection_str)
            re_field = ['re' strrep(re_str, '.', '_')]; % Handle dots if any
            d_field = ['d' strrep(num2str(deflection), '-', 'm')];
            
            data_struct.(re_field).(d_field).alpha = d(:, 1);
            data_struct.(re_field).(d_field).cl = d(:, 2);
            data_struct.(re_field).(d_field).cm = d(:, 5);
            data_struct.(re_field).(d_field).deflection = deflection;
            data_struct.(re_field).(d_field).re = str2double(strrep(re_str, 'e', 'e')); % Keep as double
            
        catch ME
            fprintf('Error reading %s: %s\n', fname, ME.message);
        end
    end

    % Get unique Reynolds numbers
    re_fields = fieldnames(data_struct);
    
    % Create one figure per Reynolds number
    for r = 1:length(re_fields)
        re_field = re_fields{r};
        d_fields = fieldnames(data_struct.(re_field));
        
        % Sort deflections for legend order
        defs = zeros(length(d_fields), 1);
        for i = 1:length(d_fields)
            defs(i) = data_struct.(re_field).(d_fields{i}).deflection;
        end
        [~, idx] = sort(defs);
        sorted_d_fields = d_fields(idx);
        
        % Create figure
        re_val = data_struct.(re_field).(d_fields{1}).re;
        fig_title = sprintf('Windex Wing Polars at Re = %.1e', re_val);
        figure('Name', fig_title, 'Color', 'w', 'Position', [100 + r*50, 100, 1200, 500]);
        
        % CL vs Alpha
        subplot(1, 2, 1);
        hold on; grid on; box on;
        xlabel('\alpha [deg]');
        ylabel('C_L');
        title('Lift Coefficient');
        
        % CM vs Alpha
        subplot(1, 2, 2);
        hold on; grid on; box on;
        xlabel('\alpha [deg]');
        ylabel('C_m');
        title('Moment Coefficient');
        
        % Plot each deflection
        colors = lines(length(sorted_d_fields));
        for i = 1:length(sorted_d_fields)
            df = sorted_d_fields{i};
            entry = data_struct.(re_field).(df);
            
            label = sprintf('d = %d^{\\circ}', entry.deflection);
            
            subplot(1, 2, 1);
            plot(entry.alpha, entry.cl, 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', label);
            
            subplot(1, 2, 2);
            plot(entry.alpha, entry.cm, 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', label);
        end
        
        subplot(1, 2, 1); legend('Location', 'best');
        subplot(1, 2, 2); legend('Location', 'best');
        sgtitle(fig_title);
    end
end
