function [Cz, Cdp, Cm, Cza] = windex_fusewingfoils(acg, fs, localpha, reynolds)
% [Cz, Cdp, Cm, Cza] = windex_foils(acg, fs, localpha, reynolds)
% 
% Determine airfoil force coefficients for each vortex segment. 
% localalpha : the "AoA" in the local coordinate system of each segment  
% reynolds   : reynolds number for each segment according to local chord and 
%              total velocity at that point 
% fs         : flight state, above all:
% fs.delta_flap : landing flap deflection (positive downward)
% fs.delta_aileron : aileron deflection, right-down is positive
% fs.delta_elevator : elevator, full HTP span 
% fs.delta_rudder : rudder, full VTP span 

  % set flap deflections from fs
  nvx = size(acg.pa,1); 
  delta = zeros(nvx,1);
  
  % shortcuts for the index sets that tell us which row belongs to which 
  % part of the aircraft 
  rw = acg.vxwing;
  
  % landing flap 

  delta(rw) = fs.delta_flap;
  
  % global data for the entire a/c
  Cz = zeros(nvx,1);
  Cdp = zeros(nvx,1);
  Cm = zeros(nvx,1);
  
  if nargout < 4
    
    % for the wing, lookup airfoil data in mwafgrid 
    [Cz(rw), Cdp(rw), Cm(rw)] = eval_afgrid(acg.mwafgrid, ...
                                localpha(rw), reynolds(rw), delta(rw));
  else
    % compute the derivatives as well
    Cza = zeros(nvx,1);
    [Cz(rw), Cdp(rw), Cm(rw), Cza(rw)] = eval_afgrid(acg.mwafgrid,...
                                 localpha(rw), reynolds(rw), delta(rw));
  end
    
end
