function [results] = generate_3d_fsim_data(V, alpha_deg)
% GENERATE_3D_FSIM_DATA
% Estimates control effectiveness derivatives (ailerons, rudder) and flap effects.
%
% Usage:
%   [results] = generate_3d_fsim_data(30, 2)
%
% V: Airspeed [m/s]
% alpha_deg: Reference Angle of Attack [deg]

    if nargin < 1; V = 30; end
    if nargin < 2; alpha_deg = 2; end

    addpath('../matlab');
    acg = make_windex();
    
    alpha_rad = alpha_deg * pi/180;
    
    % Reference Flight State (all controls zero)
    fs_ref = flight_state(acg, V, 0, alpha_rad, 0, [0 0 0]);
    fs_ref.delta_elevator = 0;
    fs_ref.delta_aileron = 0;
    fs_ref.delta_rudder = 0;
    fs_ref.delta_flap = 0;
    
    fs_ref.gamma = pointsolve(acg, fs_ref);
    [CL0, CD0, CY0, Cm0, Cl0, Cn0] = coefficients(acg, fs_ref);
    
    % Perturbation Magnitudes
    d_ctrl = 1.0 * pi/180; % 1 degree deflection
    
    % --- Aileron Derivatives ---
    % Standard convention: +delta_aileron -> Left Aileron Down, Right Aileron Up 
    % This generates a Positive Roll moment (Cl) in our windex_foils.m implementation.
    fprintf('  Estimating Aileron effectiveness...\n');
    fs_ap = fs_ref; fs_ap.delta_aileron = d_ctrl;
    fs_am = fs_ref; fs_am.delta_aileron = -d_ctrl;
    
    [CL_ap, CD_ap, CY_ap, Cm_ap, Cl_ap, Cn_ap] = run_point(acg, fs_ap);
    [CL_am, CD_am, CY_am, Cm_am, Cl_am, Cn_am] = run_point(acg, fs_am);
    
    Clda = (Cl_ap - Cl_am) / (2 * d_ctrl); % Roll due to aileron
    Cnda = (Cn_ap - Cn_am) / (2 * d_ctrl); % Yaw due to aileron (Adverse yaw)
    Cyda = (CY_ap - CY_am) / (2 * d_ctrl); % Sideforce due to aileron
    
    % --- Rudder Derivatives ---
    % Standard convention: +delta_rudder -> Trailing edge left -> Force right -> Nose right (+)
    fprintf('  Estimating Rudder effectiveness...\n');
    fs_rp = fs_ref; fs_rp.delta_rudder = d_ctrl;
    fs_rm = fs_ref; fs_rm.delta_rudder = -d_ctrl;
    
    [~, ~, CY_rp, ~, Cl_rp, Cn_rp] = run_point(acg, fs_rp);
    [~, ~, CY_rm, ~, Cl_rm, Cn_rm] = run_point(acg, fs_rm);
    
    Cndr = (Cn_rp - Cn_rm) / (2 * d_ctrl); % Yaw due to rudder
    Cldr = (Cl_rp - Cl_rm) / (2 * d_ctrl); % Roll due to rudder (Proverse/Adverse coupling)
    Cydr = (CY_rp - CY_rm) / (2 * d_ctrl); % Sideforce due to rudder
    
    % --- Flap Effects ---
    % Flaps are symmetric. Usually reported as delta_CL, delta_CD, delta_Cm per deflection.
    fprintf('  Estimating Flap effects (for 10 deg)...\n');
    df = 10 * pi/180;
    fs_f = fs_ref; fs_f.delta_flap = df;
    [CL_f, CD_f, ~, Cm_f, ~, ~] = run_point(acg, fs_f);
    
    dCL_df = (CL_f - CL0) / df;
    dCD_df = (CD_f - CD0) / df;
    dCm_df = (Cm_f - Cm0) / df;

    % Results
    results.Clda = Clda; results.Cnda = Cnda; results.Cyda = Cyda;
    results.Cndr = Cndr; results.Cldr = Cldr; results.Cydr = Cydr;
    results.dCL_df = dCL_df; results.dCD_df = dCD_df; results.dCm_df = dCm_df;

    % Output for copy-pasting
    fprintf('\n--- COPY PASTE DATA INTO make_fsim.m ---\n\n');
    fprintf('  %% Control Derivatives (per rad)\n');
    fprintf('  fsm.Clda = %.4f;  fsm.Cnda = %.4f;  fsm.Cyda = %.4f;\n', Clda, Cnda, Cyda);
    fprintf('  fsm.Cndr = %.4f;  fsm.Cldr = %.4f;  fsm.Cydr = %.4f;\n', Cndr, Cldr, Cydr);
    fprintf('\n  %% Flap Derivatives (per rad)\n');
    fprintf('  fsm.dCL_df = %.4f;  fsm.dCD_df = %.4f;  fsm.dCm_df = %.4f;\n', dCL_df, dCD_df, dCm_df);

end

function [CL, CD, CY, Cm, Cl, Cn] = run_point(acg, fs)
    fs.gamma = pointsolve(acg, fs);
    [CL, CD, CY, Cm, Cl, Cn] = coefficients(acg, fs);
end
