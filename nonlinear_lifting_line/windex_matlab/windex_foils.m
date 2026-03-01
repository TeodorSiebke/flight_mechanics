function [Cz, Cdp, Cm, Cza] = windex_foils(acg, fs, localpha, reynolds)
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
  rt = acg.vxhtp;
  rv = acg.vxvtp;
  
  ral = acg.vxail_l; % Left Aileron
  rar = acg.vxail_r; % Right Aileron
  rfl = acg.vxflap;  % Flaps
  
  % landing flap 
  delta(rfl) = fs.delta_flap;
  
  % aileron 
  % Note: check sign convention. Usually right aileron up -> positive roll?
  % fs.delta_aileron usually means "stick right". 
  % Stick right: Left aileron DOWN (+), Right aileron UP (-)? 
  % Wait, standard convention: 
  % Positive Aileron Deflection -> Right Aileron DOWN, Left Aileron UP -> Negative Roll?
  % Let's assume fs.delta_aileron is antisymmetric deflection.
  % Let's define: +delta_aileron = Left Down (+), Right Up (-) -> Roll Right? No, Left Down increases lift -> Right Roll.
  % Reference: "Right-down is positive" in comments line 10.
  % If "Right-down is positive", then for right roll (stick right), we want Left Down (+), Right Up (-).
  % Let's stick to the comment "right-down is positive" for the variable description, but for the stick input:
  % Let's assume fs.delta_aileron is the MAGNITUDE of deflection.
  % Left Aileron: +fs.delta_aileron
  % Right Aileron: -fs.delta_aileron
  
  delta( ral ) =  + fs.delta_aileron;
  delta( rar ) =  - fs.delta_aileron;

  % elevator and rudder along the span of the corresponding surfaces 
  delta( rt ) = fs.delta_elevator;
  delta( rv ) = fs.delta_rudder;
  
  % global data for the entire a/c
  Cz = zeros(nvx,1);
  Cdp = zeros(nvx,1);
  Cm = zeros(nvx,1);
  
  if nargout < 4
    
    % for the wing, lookup airfoil data in mwafgrid 
    [Cz(rw), Cdp(rw), Cm(rw)] = eval_afgrid(acg.mwafgrid, ...
                                localpha(rw), reynolds(rw), delta(rw));

    % for the tail, use htpafgrid instead - different airfoil!
    [Cz(rt), Cdp(rt), Cm(rt)] = eval_afgrid(acg.htpafgrid, ...
                                         localpha(rt), reynolds(rt), delta(rt));
    [Cz(rv), Cdp(rv), Cm(rv)] = eval_afgrid(acg.vtpafgrid, ...
                                         localpha(rv), reynolds(rv), delta(rv));
  else
    % compute the derivatives as well
    Cza = zeros(nvx,1);
    [Cz(rw), Cdp(rw), Cm(rw), Cza(rw)] = eval_afgrid(acg.mwafgrid,...
                                 localpha(rw), reynolds(rw), delta(rw));
    [Cz(rt), Cdp(rt), Cm(rt), Cza(rt)] = eval_afgrid(acg.htpafgrid, ...
                                         localpha(rt), reynolds(rt), delta(rt));
    [Cz(rv), Cdp(rv), Cm(rv), Cza(rv)] = eval_afgrid(acg.vtpafgrid, ...
                                         localpha(rv), reynolds(rv), delta(rv));
  end
    
end
