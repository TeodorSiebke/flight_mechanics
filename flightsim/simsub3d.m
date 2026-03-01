function [xdot] = simsub3d(t, x, fsm)
% simsub3d.m
% Wrapper for 3D ODE integration. Maps 12 states + 4 controls.

    % Initialize full state vector for fplmod3d (16 elements)
    xs = zeros(16,1);
    
    % Map the 12 dynamic states from ODE
    xs(1:12) = x(1:12);
    
    % Interpolate control inputs from the fsm.ctrl_set table
    % fsm.ctrl_set: [Time, de, da, dr, dp]
    ctrl = interp1(fsm.ctrl_set(:,1), fsm.ctrl_set(:,2:5), t);
    
    xs(13) = ctrl(1); % delta_e
    xs(14) = ctrl(2); % delta_a
    xs(15) = ctrl(3); % delta_r
    xs(16) = ctrl(4); % delta_p (throttle)
    
    % Optional: Add flap if needed (xs(17))
    if isfield(fsm, 'delta_flap_static')
        xs(17) = fsm.delta_flap_static;
    end

    % Calculate derivatives using the 3D model
    [xdot_full] = fplmod3d(t, xs, fsm);
    
    % Return the 12 state derivatives to the ODE solver
    xdot = xdot_full(1:12);
    xdot = xdot(:); % Enforce column vector
end
