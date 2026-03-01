%% Task 4: Trim Analysis - Linear Approximation
% Integrates Wing and Tail data to solve for the Trimmed Lift Curve using Linear Derivatives.
% Visualizes tail aerodynamics in both 2D (deflected polars) and 3D (surface).
clear; clc; close all;

% Load Configuration
if ~exist('Config', 'var')
    windex_config;
end

fprintf('\n=== TASK 4: TRIM ANALYSIS (LINEAR REVERT + 2D/3D PLOTS) ===\n');

% 1. UNPACK CONFIG
A = Config.Aero;
files = Config.Files.TrimData;
de_deg_config = Config.Files.ElevatorAngles;
de_rad_config = deg2rad(de_deg_config);

% Derived constants
Vh = (A.l_t * A.S_tail) / (A.S_wing * A.MAC); % Tail Volume Coefficient

% 2. LOAD AND STRUCTURE TAIL DATA
fprintf('Reading %d tail files...\n', length(files));

% Data containers for Linear Analysis
cl_t_at_alpha0 = zeros(size(de_deg_config));
cm_t_at_alpha0 = zeros(size(de_deg_config));
cl_alpha_t_vec = zeros(size(de_deg_config));

% Data containers for 3D Visualization
alpha_all = [];
de_all    = [];
cl_t_all  = [];

% Store for verification plotting
raw_data_alpha = cell(1, length(files));
raw_data_cl = cell(1, length(files));

for i = 1:length(files)
    currentFile = files{i};
    currentDe   = de_deg_config(i);
    
    if ~isfile(currentFile)
        warning('File not found: %s. Skipping.', currentFile);
        continue;
    end
    
    opts = detectImportOptions(currentFile);
    opts.VariableNamingRule = 'preserve';
    data = readtable(currentFile, opts);
    
    % Find Alpha Column
    varNames = data.Properties.VariableNames;
    if ismember('x____', varNames)
        alpha_loc = data.("x____");
    elseif ismember('alpha', varNames)
        alpha_loc = data.alpha;
    elseif any(contains(varNames, 'α'))
        idx = find(contains(varNames, 'α'), 1);
        alpha_loc = table2array(data(:, idx));
    elseif width(data) >= 2
        alpha_loc = table2array(data(:,2));
    else
        alpha_loc = table2array(data(:,1));
    end
    
    cl_loc = data.CL;
    cm_loc = data.Cm;
    
    % --- Accumulate for 3D Plot ---
    n_pts = length(alpha_loc);
    alpha_all = [alpha_all; alpha_loc];
    de_all    = [de_all;    repmat(currentDe, n_pts, 1)];
    cl_t_all  = [cl_t_all;  cl_loc];
    
    % --- Process for Linear Analysis ---
    % Get tail lift slope (dCl_t / d_alpha_local)
    fit_t = polyfit(deg2rad(alpha_loc), cl_loc, 1);
    cl_alpha_t_vec(i) = fit_t(1);
    
    % Get tail values at alpha_local = 0
    try
        cl_t_at_alpha0(i) = interp1(alpha_loc, cl_loc, 0, 'linear', 'extrap');
        cm_t_at_alpha0(i) = interp1(alpha_loc, cm_loc, 0, 'linear', 'extrap');
    catch
        warning('Could not interpolate alpha=0 for file %s', currentFile);
    end
    
    % Store for plotting
    raw_data_alpha{i} = alpha_loc;
    raw_data_cl{i}    = cl_loc;
end

% 3. VISUALIZATION 1: 2D Deflected Polars (Restored)
fprintf('Generating 2D Deflected Polars Plot...\n');
figure('Name', 'Deflected Elevator Polars', 'Color', 'w');
hold on;
colors = lines(length(files));
for i = 1:length(files)
    plot(raw_data_alpha{i}, raw_data_cl{i}, 'o-', 'Color', colors(i,:), ...
         'LineWidth', 1.5, 'DisplayName', sprintf('\\delta_e = %d^\\circ', de_deg_config(i)));
end
grid on; xlabel('Local Tail \alpha (deg)'); ylabel('C_{L,tail}');
title('Tail Lift Curves for Various Elevator Deflections');
legend('show', 'Location', 'best');
hold off;

% 4. VISUALIZATION 2: 3D Surface Plot (Preserved)
fprintf('Generating 3D Surface Plot...\n');
figure('Name', 'Tail Lift Data Surface', 'Color', 'w');
F_CL_tail = scatteredInterpolant(alpha_all, de_all, cl_t_all, 'linear', 'linear');
[Aq, Dq] = meshgrid(linspace(min(alpha_all), max(alpha_all), 30), ...
                    linspace(min(de_all), max(de_all), 30));
Cq = F_CL_tail(Aq, Dq);
surf(Aq, Dq, Cq, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
hold on;
plot3(alpha_all, de_all, cl_t_all, 'r.', 'MarkerSize', 10);
xlabel('Alpha (deg)'); ylabel('Elevator (deg)'); zlabel('CL_{tail}');
title('Tail Lift Characteristics Surface');
colorbar;
grid on; view(135, 30);


% 5. LINEAR DERIVATIVES CALCULATION
fprintf('Calculating Linear Derivatives...\n');

CL_alpha_t = mean(cl_alpha_t_vec);        % Tail Lift Slope (Average)
fit_cl_de  = polyfit(de_rad_config, cl_t_at_alpha0, 1);
CL_t_de    = fit_cl_de(1);                % Tail Lift change w.r.t elevator
fit_cm_de  = polyfit(de_rad_config, cm_t_at_alpha0, 1);
Cm_t_de    = fit_cm_de(1);                % Tail Moment change w.r.t elevator

% Total Aircraft Derivatives
% Total CL_alpha
CL_alpha_total = A.CL_alpha_wb + A.eta_ht * (A.S_tail/A.S_wing) * CL_alpha_t * (1 - A.deps_dalpha);

% Total Cm_alpha
Cm_alpha_total = -Vh * A.eta_ht * CL_alpha_t * (1 - A.deps_dalpha);

% Elevator Derivatives for Total Aircraft
CL_de_total = A.eta_ht * (A.S_tail/A.S_wing) * CL_t_de;
Cm_de_total = -Vh * A.eta_ht * CL_t_de; 

% Cm0 Calculation
Cm0_tail = -Vh * A.eta_ht * CL_alpha_t * (-A.eps_0); 
Cm0_total = A.Cm0_wingbody + Cm0_tail;

% 6. TRIMMED CURVE CALCULATION (Analytical)
% PHYSICS NOTE:
% The Trimmed Lift Coefficient is the sum of the Wing-Body lift and the Tail lift
% at the trim condition:
% CL_trim = CL_WB + CL_tail_trim  (Scaled by Areas and Dynamic Pressure)
%
% Usually, the tail produces a downward force (negative lift) to balance the
% natural nose-down pitching moment of the main wing (Cm0_wb < 0).
% Therefore, CL_trim is often slightly LESS than the Wing-Body lift alone.

det_val = CL_alpha_total * Cm_de_total - CL_de_total * Cm_alpha_total;

% Trimmed Lift Curve Slope
CL_alpha_trim = CL_alpha_total - (CL_de_total / Cm_de_total) * Cm_alpha_total;

% Plotting Range
a_deg_plot = linspace(Config.Sim.PlotAlphaRange(1), Config.Sim.PlotAlphaRange(2), 50);
a_rad_plot = deg2rad(a_deg_plot);

% Trimmed Lift Curve Equation (Linear)
CL_trim = -(Cm0_total * CL_de_total / Cm_de_total) + (det_val / Cm_de_total) * a_rad_plot;

% Untrimmed Curves (Fixed elevator = 0)
% Simplified linear reference for untrimmed
CL_untrimmed = CL_alpha_total * a_rad_plot + (A.eta_ht * (A.S_tail/A.S_wing) * interp1(de_deg_config, cl_t_at_alpha0, 0, 'linear', 'extrap'));

% 7. OUTPUT & PLOTTING
fprintf('--------------------------------------------------\n');
fprintf('Trimmed CL_alpha:  %.4f /rad\n', CL_alpha_trim);
fprintf('Total Cm0:         %.4f\n', Cm0_total);
fprintf('Determinant:       %.4f\n', det_val);
fprintf('--------------------------------------------------\n');

hFig = figure('Name', 'Task 4: Trim Analysis Results (Linear)', 'Color', 'w');
plot(a_deg_plot, CL_trim, 'b-', 'LineWidth', 2); hold on;
plot(a_deg_plot, CL_untrimmed, 'r--', 'LineWidth', 1.5);
plot(a_deg_plot, A.CL_alpha_wb * a_rad_plot, 'k:', 'LineWidth', 1);

grid on; xlabel('\alpha (degrees)'); ylabel('C_L');
title('Trimmed vs Untrimmed Lift Curves (Linear fit)');
legend('Trimmed Lift (C_m = 0)', 'Untrimmed (Fixed \delta_e=0)', 'Wing Body Only', 'Location', 'best');

% Save Figure
saveas(hFig, 'Trimmed_Lift_Curve_Linear.png');
fprintf('Figure saved: Trimmed_Lift_Curve_Linear.png\n');
