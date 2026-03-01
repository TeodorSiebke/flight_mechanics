function acg = make_example(nw)

  % main wing geometry 
  b = 20.0;
  c = 2.5;
  
  % htp geometry 
  ct = 0.4*c;
  bt = 0.3*b;
  xt = 8.0;
  zt = 0.6;
  
  % store reference properties
  S = c * b;
  acg.Sref = S;
  acg.bref = b;
  acg.cref = c;
  
  % aerodynamic reference point; all moments computed about this pt.
  acg.pref = [0.25*c 0 0];

  % "mesh" : number of vortex segments for main wing and tail
  % nw = 24; 
  nt = nw/2;
  
  % twist 
  troot = 3.0*pi/180;
  ttip = 1.0*pi/180;
  
  % rectangular left wing
  acg = add_wing(acg, [0.0, -0.5*b, 0.0], [0.0, 0.0, 0.0], c, c, ...
                 ttip, troot, nw/2);
             
  % rectangular right wing
  acg = add_wing(acg, [0.0, 0.0, 0.0], [0.0, 0.5*b, 0.0], c, c, ...
                 troot, ttip, nw/2);
             
  % HTP - in neutral position 
  troot = 0.0*pi/180;
  ttip  = 0.0*pi/180;           
  acg = add_wing(acg, [xt, -0.5*bt, zt], [xt, 0.5*bt, zt], ct, ct, ...
                 troot, ttip, nt);
  
  % idealized airfoil everywhere - this must be exchanged for a function 
  % that takes measured/computed nonlinear section behaviour into account!
  acg.foil = @simplefoil;
  
  % fuselage + VTP drag in the form S*CD
  acg.ffsl = 0.2;
  
end 

