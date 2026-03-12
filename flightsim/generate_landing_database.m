function [valpha, vde, cldat_land, cddat_land, cmdat_land] = generate_landing_database()
% GENERATE_LANDING_DATABASE
% Generates a 2D lookup table for CL, CD, and Cm based on Alpha and Elevator
% specifically for the LANDING configuration (20 deg flaps).

    addpath('../matlab');
    acg = make_windex();
    
    % Set ARP to 0.25 (absolute meters, matching make_fsim.m template)
    acg.pref = [0.25, 0.0, 0.0]; 

    valpha = (-5:1:20)' * pi/180; % Extended to see stall
    vde = [-5 0 5]' * pi/180;
    
    na = length(valpha);
    nd = length(vde);
    
    cldat_land = zeros(na, nd);
    cddat_land = zeros(na, nd);
    cmdat_land = zeros(na, nd);
    
    uoo = 25; % Landing speed
    fs = flight_state(acg, uoo, 0, 0, 0);
    fs.delta_flap = 20 * pi/180; % 20 deg flaps
    fs.delta_aileron = 0; fs.delta_rudder = 0;

    fprintf('Generating LANDING database (20 deg flaps) Alpha x DE [%d x %d]...\n', na, nd);
    
    for j = 1:nd
        fprintf('  Simulating Elevator: %.1f deg\n', vde(j)*180/pi);
        gamma = zeros(size(acg.pa,1), 1); % reset sweep
        for i = 1:na
            fsi = fs;
            fsi.alpha = valpha(i);
            fsi.delta_elevator = vde(j);
            fsi.vb = body_velocity(acg, uoo, fsi.alpha, 0, [0 0 0]);
            
            try
                gamma = pointsolve(acg, fsi, gamma);
                fsi.gamma = gamma;
                [cldat_land(i,j), cddat_land(i,j), ~, cmdat_land(i,j)] = coefficients(acg, fsi);
            catch
                fprintf('    Failed at Alpha=%.1f\n', valpha(i)*180/pi);
                cldat_land(i,j) = NaN;
            end
        end
    end
    
    % Format for copy-pasting
    fprintf('\n--- COPY PASTE LANDING DATA INTO make_fsim.m ---\n\n');
    
    fprintf('  fsm.cldat_land = [\n');
    for i=1:na; fprintf('    %s\n', num2str(cldat_land(i,:), '%.4f  ')); end
    fprintf('    ];\n\n');
    
    fprintf('  fsm.cddat_land = [\n');
    for i=1:na; fprintf('    %s\n', num2str(cddat_land(i,:), '%.4f  ')); end
    fprintf('    ];\n\n');
    
    fprintf('  fsm.cmdat_land = [\n');
    for i=1:na; fprintf('    %s\n', num2str(cmdat_land(i,:), '%.4f  ')); end
    fprintf('    ];\n\n');

end
