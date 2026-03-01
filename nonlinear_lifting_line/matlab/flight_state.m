function fs = flight_state(acg, uoo, alt, alpha, beta, omega)
% fs = flight_state(acg, uoo, alt, alpha, beta, omega)
% 
% acg: Initialized struct that describes the vortex geometry 
% uoo: Scalar norm of the velocity of the reference point acg.pref 
% alt: Altitude [m]
% alpha: Angle of attack [rad]
% beta: Sideslip angle [rad], zero if left out 
% omega: Vector of size (1x3) of rotation rates [rad/s]
%        omega(1) > 0 means right-wing-up roll,
%        omega(2) > 0 means pitch-up,
%        omega(3) > 0 is nose-left yaw.

    if nargin < 6
      omega = [0 0 0];
      if nargin < 5 
        beta = 0;
      end 
    end 
    
    % sanity check 
    if abs(alpha) > 30*pi/180 
      warning('Alpha looks too large.');
    end 
    if abs(beta) > 30*pi/180 
      warning('Beta looks too large.');
    end 
    if uoo == 0
      error('Sorry, refusing to fly at airspeed zero.');
    end
    
    % reference values
    fs.alpha = alpha;
    fs.beta = beta;
    fs.uref = uoo;
    fs.omega = omega;
    
    % octave-compatible replacement for atmosisa
    Tground = 15; % ground temp in Celsius
    [rho,a,T,P] = isaatm(alt, Tground-15);
    % [T, a, P, rho] = atmosisa(alt);
    
    fs.alt = alt;
    fs.rho = rho;
    fs.poo = P;
    fs.aoo = a;
    
    % compute kinematic viscosity for this altitude 
    % and ground temperatur offset Tplus
    C = 120;
    To = 291.15;
    muo = 18.27e-6;
    mu = muo * (To + C)/(T + C) * (T/To)^1.5;
    
    % dynamic viscosity 
    fs.nu = mu / rho;
    
    % velocities of the vortex segment centers 
    fs.vb = body_velocity(acg, uoo, alpha, beta, omega);
    
    % define flap deflection for every vortex element 
    fs.delta = zeros(size(acg.pa,1),1);
    
    % Default control surface deflections
    fs.delta_flap = 0;
    fs.delta_elevator = 0;
    fs.delta_aileron = 0;
    fs.delta_rudder = 0;
    
end