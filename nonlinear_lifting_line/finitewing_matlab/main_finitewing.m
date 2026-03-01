%
% Main NLLLT analysis for the Finite Wing
% Compares theoretical NLLLT results with Wind Tunnel data
%

clear all;
close all;

% --- Path Setup ---
% Add the core NLLLT library to the path
addpath('../matlab');
% Add the wind tunnel data scripts to the path
addpath('../../finite_wing/Matlab_tunneldata');

% --- Initialize Wing Model ---
fprintf('Initializing Finite Wing geometry...\n');
acg = make_finitewing();

% --- Analysis Controls ---
uoo = 30;              % Speed to analyze (m/s)
alpha_limit_min = -15;  % Minimum Alpha to analyze (deg)
alpha_limit_max = 25;  % Maximum Alpha to analyze (deg)
alpha_offset = -2.5;   % Offset to align WT data with XFOIL (deg)
alpha_sim = [alpha_limit_min:0.5:alpha_limit_max]; % Simulation alpha vector

% Check if airfoil grid is available
if isempty(acg.mwafgrid)
    error('Airfoil grid is empty. Please generate XFOIL polars.');
end

% --- Load Wind Tunnel Data ---
fprintf('Loading Wind Tunnel data (u = %d m/s)...\n', uoo);
wt_file = sprintf('../../finite_wing/Matlab_tunneldata/u%d.txt', uoo);
if exist(wt_file, 'file') ~= 2
    error('Wind tunnel data file not found: %s', wt_file);
end

[abq, bload, wload, wt_coeff, span, Sref, chord] = readwingloads(wt_file);

wt_alpha = abq(:,1) + alpha_offset;
wt_CL = wt_coeff(:,3);
wt_CD = wt_coeff(:,1);
wt_Cm = wt_coeff(:,5);

% --- Run NLLLT Solver ---
fprintf('Running NLLLT solver...');
n_pts = length(alpha_sim);
nlllt_CL = zeros(n_pts, 1);
nlllt_CD = zeros(n_pts, 1);
nlllt_Cm = zeros(n_pts, 1);

gamma_guess = zeros(size(acg.pa,1), 1);

for i = 1:n_pts
    alpha_rad = alpha_sim(i) * pi/180;
    fs = flight_state(acg, uoo, 0.0, alpha_rad, 0.0);
    
    try
        [gamma, info] = pointsolve(acg, fs, gamma_guess);
        if info <= 0, error('No convergence'); end
        fs.gamma = gamma;
        gamma_guess = fs.gamma;
        [CL, CD, CM] = postprocess(acg, fs);
        nlllt_CL(i) = CL;
        nlllt_CD(i) = CD;
        nlllt_Cm(i) = CM;
    catch
        nlllt_CL(i) = NaN;
        nlllt_CD(i) = NaN;
        nlllt_Cm(i) = NaN;
    end
end
fprintf(' Done.\n');

% --- Plotting Results ---
fprintf('Generating comparison plots...\n');

% 1. Lift Curve
figure(1);
plot(wt_alpha, wt_CL, 'ko', 'MarkerSize', 4); hold on;
plot(alpha_sim, nlllt_CL, 'r-', 'LineWidth', 1.5);
xlabel('Alpha [deg]'); ylabel('C_L'); title(sprintf('Lift Curve (u = %d m/s)', uoo));
grid on; legend('Wind Tunnel', 'NLLLT Theory', 'Location', 'NorthWest');

% 2. Drag Polar
figure(2);
plot(wt_CD, wt_CL, 'ko', 'MarkerSize', 4); hold on;
plot(nlllt_CD, nlllt_CL, 'r-', 'LineWidth', 1.5);
xlabel('C_D'); ylabel('C_L'); title(sprintf('Drag Polar (u = %d m/s)', uoo));
grid on; legend('Wind Tunnel', 'NLLLT Theory', 'Location', 'SouthEast');

% 3. Moment Coefficient
figure(3);
plot(wt_alpha, wt_Cm, 'ko', 'MarkerSize', 4); hold on;
plot(alpha_sim, nlllt_Cm, 'r-', 'LineWidth', 1.5);
xlabel('Alpha [deg]'); ylabel('C_m'); title(sprintf('Moment Curve (u = %d m/s)', uoo));
grid on; legend('Wind Tunnel', 'NLLLT Theory', 'Location', 'SouthWest');

fprintf('Analysis complete.\n');
