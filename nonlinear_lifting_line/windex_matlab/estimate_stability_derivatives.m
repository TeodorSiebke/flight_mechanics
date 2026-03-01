function [results, table_out] = estimate_stability_derivatives(V, alpha_deg)
% ESTIMATE_STABILITY_DERIVATIVES
% Estimates stability derivatives (beta, p, q, r) for the Windex aircraft.
%
% Usage:
%   [results, table] = estimate_stability_derivatives(30, 2)
%
% V: Airspeed [m/s]
% alpha_deg: Reference Angle of Attack [deg]

    if nargin < 1; V = 30; end
    if nargin < 2; alpha_deg = 2; end

    addpath('../matlab');
    acg = make_windex();
    
    alpha_rad = alpha_deg * pi/180;
    
    % Reference Flight State
    fs_ref = flight_state(acg, V, 0, alpha_rad, 0, [0 0 0]);
    fs_ref.gamma = pointsolve(acg, fs_ref);
    [CL0, CD0, CY0, Cm0, Cl0, Cn0] = coefficients(acg, fs_ref);
    
    % Perturbation Magnitudes (Reduced for better linearity)
    d_beta = 0.1 * pi/180;    % 0.1 degree
    d_p = 0.01;               % rad/s
    d_q = 0.01;               % rad/s
    d_r = 0.01;               % rad/s
    
    % Non-dimensionalizing Factors
    b = acg.bref;
    c = acg.cref;
    p_hat_factor = b / (2*V);
    q_hat_factor = c / (2*V);
    r_hat_factor = b / (2*V);

    % --- Sideslip Derivatives (Beta) ---
    fprintf('  Calculating Beta derivatives...\n');
    [cl_bp, cn_bp, cy_bp] = run_point(acg, V, alpha_rad,  d_beta, [0 0 0], fs_ref.gamma);
    [cl_bm, cn_bm, cy_bm] = run_point(acg, V, alpha_rad, -d_beta, [0 0 0], fs_ref.gamma);
    Clb = (cl_bp - cl_bm) / (2 * d_beta);
    Cnb = (cn_bp - cn_bm) / (2 * d_beta);
    Cyb = (cy_bp - cy_bm) / (2 * d_beta);

    % --- Roll Rate Derivatives (p) ---
    fprintf('  Calculating Roll (p) derivatives...\n');
    % NLL mapping: omega(1) = -p
    [cl_pp, cn_pp, cy_pp] = run_point(acg, V, alpha_rad, 0, [-d_p 0 0], fs_ref.gamma);
    [cl_pm, cn_pm, cy_pm] = run_point(acg, V, alpha_rad, 0, [+d_p 0 0], fs_ref.gamma);
    Clp = (cl_pp - cl_pm) / (2 * d_p * p_hat_factor);
    Cnp = (cn_pp - cn_pm) / (2 * d_p * p_hat_factor);
    Cyp = (cy_pp - cy_pm) / (2 * d_p * p_hat_factor);

    % --- Pitch Rate Derivatives (q) ---
    fprintf('  Calculating Pitch (q) derivatives...\n');
    % NLL mapping: omega(2) = q
    [cl_qp, cn_qp, cy_qp, CL_qp, CD_qp, Cm_qp] = run_point(acg, V, alpha_rad, 0, [0 +d_q 0], fs_ref.gamma);
    [cl_qm, cn_qm, cy_qm, CL_qm, CD_qm, Cm_qm] = run_point(acg, V, alpha_rad, 0, [0 -d_q 0], fs_ref.gamma);
    Cmq = (Cm_qp - Cm_qm) / (2 * d_q * q_hat_factor);
    CLq = (CL_qp - CL_qm) / (2 * d_q * q_hat_factor);
    CDq = (CD_qp - CD_qm) / (2 * d_q * q_hat_factor);

    % --- Yaw Rate Derivatives (r) ---
    fprintf('  Calculating Yaw (r) derivatives...\n');
    % NLL mapping: omega(3) = -r
    [cl_rp, cn_rp, cy_rp] = run_point(acg, V, alpha_rad, 0, [0 0 -d_r], fs_ref.gamma);
    [cl_rm, cn_rm, cy_rm] = run_point(acg, V, alpha_rad, 0, [0 0 +d_r], fs_ref.gamma);
    Clr = (cl_rp - cl_rm) / (2 * d_r * r_hat_factor);
    Cnr = (cn_rp - cn_rm) / (2 * d_r * r_hat_factor);
    Cyr = (cy_rp - cy_rm) / (2 * d_r * r_hat_factor);

    % Results formatting
    results.Clb = Clb; results.Cnb = Cnb; results.Cyb = Cyb;
    results.Clp = Clp; results.Cnp = Cnp; results.Cyp = Cyp;
    results.Cmq = Cmq; results.CLq = CLq; results.CDq = CDq;
    results.Clr = Clr; results.Cnr = Cnr; results.Cyr = Cyr;

    if nargout > 1
        table_out = results; % Simple for now, can format as a table if needed
    end

    % Display Output
    fprintf('=== Windex Stability Derivatives (V=%.1f m/s, alpha=%.1f deg) ===\n', V, alpha_deg);
    fprintf('\nLongitudinal Derivatives:\n');
    fprintf('  Cmq: %8.4f\n', Cmq);
    fprintf('  CLq: %8.4f\n', CLq);
    
    fprintf('\nLateral-Directional Derivatives (per rad):\n');
    fprintf('        |    beta    |     p      |     r      \n');
    fprintf('  ---------------------------------------------\n');
    fprintf('  Cl    | %10.4f | %10.4f | %10.4f \n', Clb, Clp, Clr);
    fprintf('  Cn    | %10.4f | %10.4f | %10.4f \n', Cnb, Cnp, Cnr);
    fprintf('  CY    | %10.4f | %10.4f | %10.4f \n', Cyb, Cyp, Cyr);
    
end

function [Cl, Cn, CY, CL, CD, Cm] = run_point(acg, V, alpha, beta, omega, gamma0)
    fs = flight_state(acg, V, 0, alpha, beta, omega);
    fs.gamma = pointsolve(acg, fs, gamma0);
    [CL, CD, CY, Cm, Cl, Cn] = coefficients(acg, fs);
end
