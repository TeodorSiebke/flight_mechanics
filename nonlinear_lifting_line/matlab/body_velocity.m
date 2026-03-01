function vb = body_velocity(acg, uoo, alpha, beta, omega)
  
  
  % body velocity, speed of each vortex segment in 
  % its own coordinate system
  T = body2wind(alpha, beta);
  vxy = (T' * [uoo; 0; 0])';
  
  % contributions of the rotation rates - rotation is taken about the 
  % aerodynamic reference point set when creating acg  
  ns = size(acg.pa, 1);
  r = 0.5*(acg.pa + acg.pb) - repmat(acg.pref, [ns 1]);
% BUG fixed 20180301
%  vb = repmat(-vxy, [ns 1]); + cross(repmat(omega, [ns 1]), r);
  vb = repmat(-vxy, [ns 1]) + cross(repmat(omega, [ns 1]), r);
  
end 
