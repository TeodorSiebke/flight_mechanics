%% Windex Analysis Configuration File
% This file allows you to control all variables from one place.
% Run this script to load variables into your workspace, or let the task scripts call it.

%% 1. Geometry Definitions
% Using the format: [Y-pos, Chord, xOffset_from_RootLE]
% xOffset is relative to that surface's own Root Leading Edge

% Main Wing Geometry
Config.Geom.Wing = [
    0.000, 0.750, 0.000;   % Root
    3.500, 0.600, 0.070;   % Mid
    6.000, 0.330, 0.200    % Tip
];

% Horizontal Tail Geometry
Config.Geom.Tail = [
    0.000, 0.570, 0.000;   % Root
    1.130, 0.330, 0.200    % Tip
];

%% 2. Aircraft & Aerodynamic Constants
% Geometric reference values
Config.Aero.S_wing = 7.41;          % Wing reference area (m^2)
Config.Aero.S_tail = 0.87;          % Horizontal tail area (m^2)

% Main Wing MAC properties
Config.Aero.MAC_wing    = 0.6119;         % Main Wing MAC (m)
Config.Aero.X_LE_wing   = 0.0650;         % Main Wing MAC Leading Edge X-pos (m)

% Horizontal Tail MAC properties
Config.Aero.MAC_tail    = 0.4607;         % Horizontal Tail MAC (m)
Config.Aero.X_LE_tail   = 0.0911;         % Horizontal Tail MAC Leading Edge X-pos (m)

% Legacy support (if needed by other scripts)
Config.Aero.MAC    = Config.Aero.MAC_wing;         
Config.Aero.l_t    = 2.681;         % Distance from CG/Wing AC to Tail AC (m)

% Aerodynamic Coefficients
Config.Aero.h_nwb       = 0.2486;            % Neutral point of wing-body (fraction of MAC)
Config.Aero.CL_alpha_wb = 0.1003 * 180/pi;   % Wing-Body Lift slope (per rad)
Config.Aero.CL_alpha_ht = 0.0774 * 180/pi;   % Tail Lift slope (per rad)
Config.Aero.eta_ht      = 0.95;              % Tail efficiency factor
Config.Aero.eps_0       = 0.02;              % Downwash at zero alpha (rad)
Config.Aero.deps_dalpha = 0.176;             % Downwash gradient (d_epsilon/d_alpha)

% Missing value from original script (Placeholder)
Config.Aero.Cm0_wingbody = -0.10;            % Zero-lift moment coefficient of Wing-Body (Estimate)

%% 3. File Paths
% Define file names here so you don't have to select them every time.

% For Task 3: Polar Analysis (Linear Region)
Config.Files.PolarAnalysisWing = 'T1-20_0 m_s-LLT.csv'; 
Config.Files.PolarAnalysisTail = 'T1-20_0 m_s-LLT_HT.csv';

% For Task 4: Trim Analysis - List of Tail LLT CSVs for different elevator deflections
Config.Files.TrimData = {
    'T1-30_0 m_s-LLT_-10.csv', ...
    'T1-30_0 m_s-LLT_-5.csv',  ...
    'T1-30_0 m_s-LLT_0.csv',   ...
    'T1-30_0 m_s-LLT_5.csv',   ...
    'T1-30_0 m_s-LLT_10.csv'
};

% Corresponding elevator angles [deg]
Config.Files.ElevatorAngles = [-10, -5, 0, 5, 10];

%% 4. Simulation Settings
Config.Sim.LinearAlphaRange = [-5, 10];      % Range of alpha for linear analysis [deg]
Config.Sim.PlotAlphaRange   = [-5, 12];      % Range of alpha for plotting trim curves [deg]

disp('Windex Configuration Loaded.');
