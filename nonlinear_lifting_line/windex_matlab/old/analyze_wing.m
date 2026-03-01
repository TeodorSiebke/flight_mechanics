% analyze_wing.m
% Isolated analysis of the Windex Main Wing.
% Calculates MAC, XLE_MAC, and Neutral Point referenced to MAC LE.

function analyze_wing()
    clear all; close all;
    addpath('../matlab');

    fprintf('--- ANALYZING WINDEX MAIN WING ---\n');
    acg_full = make_windex();
    rw = acg_full.vxwing;

    % 1. Numerical Geometry Properties
    pa = acg_full.pa(rw,:);
    pb = acg_full.pb(rw,:);
    chord = acg_full.chord(rw);
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
    acg.pa = pa; acg.pb = pb;
    acg.xh = acg_full.xh(rw,:); acg.zup = acg_full.zup(rw,:);
    acg.chord = chord; acg.incidence = acg_full.incidence(rw);
    acg.vxwing = 1:numel(rw); acg.vxhtp = []; acg.vxvtp = [];
    
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
    % Cm = CL * (h_ref - h_ac) -> dCm/dCL = -h_ac (if h_ref = 0 at MAC LE)
    dcm_dcl = cma / cla;
    np_perc = -dcm_dcl * 100;

    fprintf('\nAnalysis Results (Alpha -2 to 5 deg):\n');
    fprintf('CL_alpha:        %.4f /rad\n', cla);
    fprintf('CM_alpha:        %.4f /rad\n', cma);
    fprintf('dCm/dCL:         %.4f\n', dcm_dcl);
    fprintf('Neutral Point:   %.2f%% of MAC (aft of MAC LE)\n', np_perc);

    % 5. Plotting
    figure('Name', 'Main Wing Analysis');
    subplot(2,1,1);
    plot(alpha_sweep_deg, CL, 'b-o', 'LineWidth', 1.5); hold on;
    plot(alpha_sweep_deg(lin_idx), polyval(p_cl, al), 'r--', 'LineWidth', 1);
    xlabel('\alpha [deg]'); ylabel('C_L'); grid on; title('Wing Lift Curve');
    
    subplot(2,1,2);
    plot(alpha_sweep_deg, CM, 'b-o', 'LineWidth', 1.5); hold on;
    plot(alpha_sweep_deg(lin_idx), polyval(p_cm, al), 'r--', 'LineWidth', 1);
    xlabel('\alpha [deg]'); ylabel('C_m (ref @ MAC LE)'); grid on; title('Wing Moment Curve');

end
