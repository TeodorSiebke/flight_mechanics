function [Cz, Cdp, Cm, Cza] = finitewing_foils(acg, fs, localpha, reynolds)
% [Cz, Cdp, Cm, Cza] = finitewing_foils(acg, fs, localpha, reynolds)
% 
% Determine airfoil force coefficients for each vortex segment. 

  nvx = size(acg.pa,1); 
  rw = acg.vxwing;
  
  Cz = zeros(nvx,1);
  Cdp = zeros(nvx,1);
  Cm = zeros(nvx,1);
  
  % Assuming no flaps for the finite wing windtunnel model
  delta = zeros(nvx,1);

  if isempty(acg.mwafgrid)
      error('Airfoil grid (mwafgrid) is empty. Please generate XFOIL polars.');
  end

  if nargout < 4
    [Cz(rw), Cdp(rw), Cm(rw)] = eval_afgrid(acg.mwafgrid, ...
                                localpha(rw), reynolds(rw), delta(rw));
  else
    Cza = zeros(nvx,1);
    [Cz(rw), Cdp(rw), Cm(rw), Cza(rw)] = eval_afgrid(acg.mwafgrid,...
                                 localpha(rw), reynolds(rw), delta(rw));
  end
    
end
