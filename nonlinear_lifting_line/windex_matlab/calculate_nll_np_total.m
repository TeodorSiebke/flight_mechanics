function [h_n, X_LE_wing_mac, x_np] = calculate_nll_np_total(acg)
% [h_n, X_LE_wing_mac, x_np] = calculate_nll_np_total(acg)
% Calculates the Total Aircraft Neutral Point directly from NLL simulation forces.

if nargin < 1
    acg = make_windex();
    interactive = true;
else
    interactive = false;
end

if interactive
    fprintf('=== WINDEX TOTAL NEUTRAL POINT CALCULATION (DIRECT NLL) ===\n\n');
end

% 1. SETUP MODEL
rw = acg.vxwing;
l_seg = sqrt(sum((acg.pb - acg.pa).^2, 2));
S_segs = acg.chord .* l_seg;
S_wing = sum(S_segs(rw)); 
MAC_wing = sum(acg.chord(rw).^2 .* l_seg(rw)) / S_wing; 
X_LE_wing_mac = sum(acg.pa(rw,1) .* S_segs(rw)) / S_wing;
S_ref = acg.Sref;
c_ref = acg.cref; 

if interactive
    fprintf('Reference Geometry:\n');
    fprintf('  S_ref: %.4f m^2, c_ref: %.4f m, MAC LE: %.4f m\n', S_ref, c_ref, X_LE_wing_mac);
end

x_ref = X_LE_wing_mac;
acg.pref = [x_ref, 0, 0]; 

% 2. RUN ANALYSIS
uoo = 30; 
alpha_sweep_deg = -2:1:8; 
alpha_sweep_rad = alpha_sweep_deg * pi/180;
fs = flight_state(acg, uoo, 0, 0, 0);
fs.delta_elevator = 0; fs.delta_rudder = 0; fs.delta_flap = 0; fs.delta_aileron = 0;

[CL, CD, CM] = plainpolar(acg, fs, alpha_sweep_rad);

% 3. CALCULATE STABILITY
p_cl = polyfit(alpha_sweep_rad, CL, 1);
cla = p_cl(1); 
p_cm = polyfit(alpha_sweep_rad, CM, 1);
cma = p_cm(1); 
dcm_dcl = cma / cla;
x_np = x_ref - dcm_dcl * c_ref;
h_n = (x_np - X_LE_wing_mac) / c_ref;

if interactive
    fprintf('--- RESULTS ---\n');
    fprintf('dCL/da: %.4f, dCm/da: %.4f, dCm/dCL: %.4f\n', cla, cma, dcm_dcl);
    fprintf('X_NP: %.4f m, NP (%% MAC): %.2f%%\n', x_np, h_n * 100);
    
    figure('Name', 'Total Aircraft Stability');
    subplot(2,1,1); plot(alpha_sweep_deg, CL, 'b-o'); grid on; title('CL vs alpha');
    subplot(2,1,2); plot(CL, CM, 'r-o'); grid on; title('Cm vs CL');
end

end
