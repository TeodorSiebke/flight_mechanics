function acg = make_wtm2(nw)
% acg = make_wtm2()
%
% Create the geometry structure for the finite wing wind tunnel lab
% This function needs to be completed first!

  % main wing geometry 
  b = 
  ctip =  
  croot = 
  dbody = 
  
  % store reference properties
  S = 
  acg.Sref = S;
  acg.bref = b;
  acg.cref = 0.5*(ctip + croot);
  
  % aerodynamic reference point: use the balance center for this
  xref = 
  acg.pref = [ xref, 0.0, 0.0 ];

            
  % setting angle, in radian
  aset = 
  
  % trapezoidal left wing
  nside = nw/2 - 2;
  xle_tip = 
  xle_root = 
  acg = add_wing(acg, [xle_tip, -0.5*b, 0.0], ...
                      [xle_root, -0.5*dbody, 0.0], ctip, croot, ...
                      aset, aset, nside);
                  
  % center segment
  acg = add_wing(acg, [xle_root, -0.5*dbody, 0.0], ...
                      [xle_root, +0.5*dbody, 0.0], croot, croot, ...
                      aset, aset, 4);
                  
  % trapezoidal right wing
  acg = add_wing(acg, [xle_root, 0.5*dbody, 0.0], ...
                      [xle_tip, 0.5*b, 0.0], croot, ctip, ...
                      aset, aset, nside);
  
  % generate lookup table for H15 airfoil, using XFoil data
  acg.afgrid = build_h15();
  
  % use the lookup table 
  acg.foil = @lookupfoil;
  
  % drag contribution of the central body - better to add that later 
  % from measurements taken on the isolated body
  acg.ffsl = 0.0;
  
end 

