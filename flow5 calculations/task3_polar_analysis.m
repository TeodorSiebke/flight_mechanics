%% Task 3: Polar Analysis (Linear Region)
% Analyzes Flow5 export files for both Main Wing and Horizontal Tail.
% Calculates CL_alpha and Neutral Point (NP) location for each.
clear; clc; close all;

% Load Configuration
if ~exist('Config', 'var')
    windex_config;
end

fprintf('\n=== TASK 3: POLAR ANALYSIS ===\n');

%% 1. ANALYZE MAIN WING
fprintf('\n--- Processing Main Wing ---\n');
analyze_surface(Config.Files.PolarAnalysisWing, ...
                Config.Aero.MAC_wing, ...
                Config.Aero.X_LE_wing, ...
                'Main Wing', ...
                Config);

%% 2. ANALYZE HORIZONTAL TAIL
fprintf('\n--- Processing Horizontal Tail ---\n');
analyze_surface(Config.Files.PolarAnalysisTail, ...
                Config.Aero.MAC_tail, ...
                Config.Aero.X_LE_tail, ...
                'Horizontal Tail', ...
                Config);

%% HELPER FUNCTION
function analyze_surface(filename, MAC, X_LE_MAC, TitleString, Config)
    fprintf('Target File: %s\n', filename);
    
    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    % --- READ DATA ---
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    data = readtable(filename, opts);
    
    % Robust Column Extraction (copied from Task 4 fix)
    varNames = data.Properties.VariableNames;
    
    % 1. Extract Alpha
    if ismember('x____', varNames)
        alpha_full = data.("x____");
    elseif ismember('alpha', varNames)
        alpha_full = data.alpha;
    elseif any(contains(varNames, 'α'))
        idx = find(contains(varNames, 'α'), 1);
        alpha_full = table2array(data(:, idx));
    elseif width(data) >= 2
        alpha_full = table2array(data(:,2)); % Fallback
    else
        alpha_full = table2array(data(:,1)); % Last resort
    end
    
    % 2. Extract CL
    if ismember('CL', varNames)
        CL_full = data.CL;
    elseif ismember('Cl', varNames)
        CL_full = data.Cl;
    else
        % Fallback: usually column 5 for Flow5 LLT csv
        CL_full = table2array(data(:,5)); 
    end
    
    % 3. Extract Cm
    if ismember('Cm', varNames)
        Cm_full = data.Cm;
    elseif ismember('CMy', varNames)
        Cm_full = data.CMy;
    else
        % Fallback: usually column 10 for Flow5 LLT csv (double check?)
        % In the viewed file, Cm is column 10
        Cm_full = table2array(data(:,10)); 
    end

    % --- FILTER DATA ---
    valid_idx = ~isnan(CL_full) & ~isnan(Cm_full) & ~isnan(alpha_full);
    alpha_full = alpha_full(valid_idx);
    CL_full    = CL_full(valid_idx);
    Cm_full    = Cm_full(valid_idx);
    
    % Filter for linear region
    range = Config.Sim.LinearAlphaRange;
    linear_mask = (alpha_full >= range(1)) & (alpha_full <= range(2));
    alpha_lin = alpha_full(linear_mask);
    CL_lin    = CL_full(linear_mask);
    Cm_lin    = Cm_full(linear_mask);
    
    if isempty(alpha_lin)
        warning('No data points found in linear range [%g, %g]. Using all data.', range(1), range(2));
        alpha_lin = alpha_full; CL_lin = CL_full; Cm_lin = Cm_full;
    end
    
    % --- CALCULATIONS ---
    % A. Calculate CL_alpha
    p_lift = polyfit(deg2rad(alpha_lin), CL_lin, 1);
    CL_alpha_rad = p_lift(1);
    CL_alpha_deg = CL_alpha_rad * (pi/180);
    
    % B. Calculate Neutral Point (Aerodynamic Center)
    % Slope: dCm / dCL
    p_stab = polyfit(CL_lin, Cm_lin, 1);
    dCmdCL = p_stab(1); 
    
    % X_NP calculations (Standard formula for stability: dCm/dCL = -(X_np - X_ref)/MAC )
    % Here the reference is likely the root LE or wherever the moment was taken about.
    % Flow5 moments are usually about the reference point defined in the analysis.
    % Assuming the Moment in file is about the Root Leading Edge (default if not specified otherwise usually):
    % X_NP - X_ref = -dCmdCL * MAC
    % If X_ref = 0 (Root LE), then X_NP = -dCmdCL * MAC
    
    X_NP_from_MomentRef = -dCmdCL * MAC; 
    
    % Note: Verify if the CSV moment is about CoG or Root LE. 
    % Usually strictly aerodynamic output (Cm) is about the defined CoG in Flow5 or 0,0,0.
    % We assume here it is about the "Reference Point" which is often 0 for simple analysis
    % But let's stick to the user's previous logic: X_NP_global = -dCmdCL * MAC
    
    X_NP_global = X_NP_from_MomentRef; 
    X_NP_relative = X_NP_global - X_LE_MAC;
    X_NP_percent  = (X_NP_relative / MAC) * 100;
    
    % --- OUTPUT ---
    fprintf('  MAC used:                  %.4f m\n', MAC);
    fprintf('  X_LE of MAC:               %.4f m\n', X_LE_MAC);
    fprintf('  Lift Curve Slope (CL_a):   %.4f /deg\n', CL_alpha_deg);
    fprintf('  Stability Slope (dCm/dCL): %.4f\n', dCmdCL);
    fprintf('  NP location (from Ref):    %.4f m\n', X_NP_global);
    fprintf('  NP Position (%% MAC):       %.2f %%\n', X_NP_percent);
    
    % --- PLOTTING ---
    figure('Name', ['Task 3: ' TitleString], 'Color', 'w');
    
    % Plot 1: Lift Curve
    subplot(1,2,1); hold on; grid on;
    plot(alpha_full, CL_full, '.', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw Data');
    plot(alpha_lin, CL_lin, 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'Linear Range');
    % Plot Fit Line
    fplot(@(x) p_lift(1)*deg2rad(x) + p_lift(2), [min(alpha_full) max(alpha_full)], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');
    title([TitleString ' - Lift Curve']);
    xlabel('\alpha (deg)'); ylabel('C_L');
    legend('Location', 'best');
    
    % Plot 2: Stability Curve
    subplot(1,2,2); hold on; grid on;
    plot(CL_full, Cm_full, '.', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw Data');
    plot(CL_lin, Cm_lin, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Linear Range');
    % Plot Fit Line
    fplot(@(x) p_stab(1)*x + p_stab(2), [min(CL_full) max(CL_full)], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');
    title([TitleString ' - Stability Curve']);
    xlabel('C_L'); ylabel('C_m');
    legend('Location', 'best');
end
