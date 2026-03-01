function [gamma, info] = pointsolve(acg, fs, gamma0)
% gamma = pointsolve(acg, fs)
%
% Solve for the vortex strength gamma at the single operating point 
% specified by fs.
%
% acg: struct containing vortex geometry and airfoil data lookup function
% fs: struct describing the flight condition 
% gamma0: initial guess for the solution (can be left out) 

  if nargin < 3
    if isfield(fs, 'gamma')
      gamma0 = fs.gamma;
    else 
      gamma0 = zeros(size(acg.pa,1),1);
    end 
  end 

  fobj = @(y) residual(acg, fs, y);
  opt = optimset('Jacobian', 'on', 'Display', 'off');
  [gamma, r, info] = fsolve(fobj, gamma0, opt);
  
  % complain when not converged
  if info == 0
    error('fsolve did not converge: iteration limit exceeded.');
  elseif info == -3
    fprintf(1, 'Failed for Alpha: %.1f deg, Beta: %.1f deg.\n', ...
               fs.alpha*180/pi, fs.beta*180/pi);
    error('fsolve did not converge: trust region too small.');
  end 
  
end 