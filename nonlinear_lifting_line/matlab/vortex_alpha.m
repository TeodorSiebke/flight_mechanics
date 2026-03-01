function [alpha, vax, vup] = vortex_alpha(acg, v)
% alpha = vortex_alpha(acg, v)
%
% Determine local angle of attack as seen by each vortex segment, where AoA
% is meant in the airfoil coordinate system (so it would be "beta" for a VTP
% vortex segment). 
  
  % upward velocity component 
  vup = sum(acg.zup .* v, 2);
  
  % longitudinal velocity component 
  vax = sum(acg.xh .* v, 2);
  alpha = acg.incidence + atan2(vup, vax);
  
end