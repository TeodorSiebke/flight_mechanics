function [CL, CD, CM] = plainpolar(acg, fs, valpha)
% [CL, CD, CM] = plainpolar(acg, fs, valpha)
%
% Compute a simple polar diagram without adjusting for trim, i.e. points on
% this polar will not correspond to trimmed flight conditions.
%
% acg : aircraft geometry 
% fs  : flight condition 
% valpha : an array of AoA values at which the polar will be evaluated
%
% CL, CD, CM : Lift, drag and moment coefficients at valpha.

  na = numel(valpha);
  CL = zeros(na, 1);
  CD = zeros(na, 1);
  CM = zeros(na, 1);
  
  if isfield(fs, 'gamma')
    gamma = fs.gamma;
  else 
    gamma = zeros( size(acg.pa,1), 1 ); 
  end   
  
  fsi = fs;
  for ki = 1:na
    fsi.alpha = valpha(ki);
    fsi.vb = body_velocity(acg, fs.uref, valpha(ki), fs.beta, fs.omega);
    gamma = pointsolve(acg, fsi, gamma);
    fsi.gamma = gamma;
    [CL(ki), CD(ki), CC, CM(ki), Cl, Cn] = coefficients(acg, fsi);
  end 

end 