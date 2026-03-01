function [valpha, vde, cldat, cddat, cmdat] = generate_fsim_database(uoo, delta_flap_deg)
% GENERATE_FSIM_DATABASE
% Generates a 2D lookup table for CL, CD, and Cm based on Alpha and Elevator.

    if nargin < 1
        uoo = 30; % Default speed
    end
    if nargin < 2
        delta_flap_deg = 0;
    end

    addpath('../matlab');
    acg = make_windex();
    
    % Set ARP to 0.25 (absolute meters, matching make_fsim.m template)
    acg.pref = [0.25, 0.0, 0.0]; 

    valpha = (-5:1:15)' * pi/180;
    vde = [-5 0 5]' * pi/180;
    
    na = length(valpha);
    nd = length(vde);
    
    cldat = zeros(na, nd);
    cddat = zeros(na, nd);
    cmdat = zeros(na, nd);
    
    % uoo is now passed as an argument
    fs = flight_state(acg, uoo, 0, 0, 0);
    fs.delta_flap = delta_flap_deg * pi/180; 
    fs.delta_aileron = 0; fs.delta_rudder = 0;

    fprintf('Generating database Alpha x DE [%d x %d] for delta_flap = %.1f deg...\n', na, nd, delta_flap_deg);
    
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
                [cldat(i,j), cddat(i,j), ~, cmdat(i,j)] = coefficients(acg, fsi);
            catch
                fprintf('    Failed at Alpha=%.1f\n', valpha(i)*180/pi);
                cldat(i,j) = NaN;
            end
        end
    end
    
    % Format for copy-pasting
    fprintf('\n--- COPY PASTE DATA INTO make_fsim.m ---\n\n');
    
    fprintf('  fsm.valpha = [\n');
    fprintf('    %.4f\n', valpha);
    fprintf('    ];\n\n');
    
    fprintf('  fsm.vde = [\n');
    fprintf('    %.4f\n', vde);
    fprintf('    ];\n\n');
    
    fprintf('  fsm.cldat = [\n');
    for i=1:na; fprintf('    %s\n', num2str(cldat(i,:), '%.4f  ')); end
    fprintf('    ];\n\n');
    
    fprintf('  fsm.cddat = [\n');
    for i=1:na; fprintf('    %s\n', num2str(cddat(i,:), '%.4f  ')); end
    fprintf('    ];\n\n');
    
    fprintf('  fsm.cmdat = [\n');
    for i=1:na; fprintf('    %s\n', num2str(cmdat(i,:), '%.4f  ')); end
    fprintf('    ];\n\n');

end
