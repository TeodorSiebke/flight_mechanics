% analyze_tail.m
% Isolated analysis of the Windex Horizontal Tail.
% Calculates MAC, XLE_MAC, and Neutral Point referenced to MAC LE.

function analyze_tail()
    clear all; close all;
    addpath('../matlab');

    fprintf('--- ANALYZING WINDEX HORIZONTAL TAIL ---\n');
    acg_full = make_windex();
    rt = acg_full.vxhtp;

    % 1. Numerical Geometry Properties
    pa = acg_full.pa(rt,:);
    pb = acg_full.pb(rt,:);
    chord = acg_full.chord(rt);
    dy = abs(pb(:,2) - pa(:,2));
    
    S = sum(chord .* dy);
    mac = sum(chord.^2 .* dy) / S;
    
    % Vortex is at 25% chord, so LE = 0.5*(pa+pb) - 0.25*chord
    % We also calculate Mean Z to avoid Z-couple bias in AC detection
    p_mid = 0.5 * (pa + pb);
    xle_mac = sum((p_mid(:,1) - 0.25*chord) .* (chord .* dy)) / S;
    z_mac   = sum(p_mid(:,3) .* (chord .* dy)) / S;

    fprintf('Ref Area (S):    %.4f m^2\n', S);
    fprintf('MAC:             %.4f m\n', mac);
    fprintf('X_LE_MAC:        %.4f m\n', xle_mac);
    fprintf('Z_Ref (Mean):    %.4f m\n', z_mac);

    % 2. Setup NLLT Model
    acg = acg_full;
    acg.mwafgrid = acg_full.htpafgrid; % Use tail airfoil
    acg.pa = pa; acg.pb = pb;
    acg.xh = acg_full.xh(rt,:); acg.zup = acg_full.zup(rt,:);
    acg.chord = chord; acg.incidence = acg_full.incidence(rt);
    acg.vxwing = []; acg.vxhtp = 1:numel(rt); acg.vxvtp = [];
    
    acg.Sref = S;
    acg.cref = mac;
    acg.pref = [xle_mac, 0, z_mac]; % REFERENCE AT MAC LE, Z-ALIGNED

    % 3. Run Alpha Sweep
    uoo = 30;
    alpha_sweep_deg = -2:1:10;
    alpha_sweep_rad = alpha_sweep_deg * pi/180;
    
    fs = flight_state(acg, uoo, 0, 0, 0);
    fs.delta_elevator = 0; fs.delta_rudder = 0; fs.delta_flap = 0; fs.delta_aileron = 0;

    [CL, CD, CM] = plainpolar(acg, fs, alpha_sweep_rad);

    % 4. Linear Slope Analysis (-2 to 5 degrees)
    lin_idx = (alpha_sweep_deg >= -2) & (alpha_sweep_deg <= 5);
    al = alpha_sweep_rad(lin_idx);
    
    p_cl = polyfit(al, CL(lin_idx), 1);
    cla = p_cl(1); % dCL/da
    
    p_cm = polyfit(al, CM(lin_idx), 1);
    cma = p_cm(1); % dCM/da
    
    % Neutral Point calculation
    dcm_dcl = cma / cla;
    np_perc = -dcm_dcl * 100;

    fprintf('\nAnalysis Results (Alpha -2 to 5 deg):\n');
    fprintf('CL_alpha:        %.4f /rad\n', cla);
    fprintf('CM_alpha:        %.4f /rad\n', cma);
    fprintf('dCm/dCL:         %.4f\n', dcm_dcl);
    fprintf('Neutral Point:   %.2f%% of MAC (aft of MAC LE)\n', np_perc);

    % 5. Plotting
    figure('Name', 'Horizontal Tail Analysis');
    subplot(2,1,1);
    plot(alpha_sweep_deg, CL, 'r-s', 'LineWidth', 1.5); hold on;
    plot(alpha_sweep_deg(lin_idx), polyval(p_cl, al), 'k--', 'LineWidth', 1);
    xlabel('\alpha [deg]'); ylabel('C_L'); grid on; title('Tail Lift Curve');
    
    subplot(2,1,2);
    plot(alpha_sweep_deg, CM, 'r-s', 'LineWidth', 1.5); hold on;
    plot(alpha_sweep_deg(lin_idx), polyval(p_cm, al), 'k--', 'LineWidth', 1);
    xlabel('\alpha [deg]'); ylabel('C_m (ref @ MAC LE)'); grid on; title('Tail Moment Curve');

end
