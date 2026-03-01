function acg = make_finitewing()
% store reference properties for the Finite Wing
% Data from finite_wing/liftline.m

  b = 1.125;     % Total span (m)
  bw = 0.125;    % Width of constant chord mid section (m)
  croot = 0.277; % Root section chord (m)
  ctip = 0.21;   % Tip chord (m)
  
  % Reference area (trapezoidal + rectangular)
  Sref = (croot + ctip) * (b - bw) / 2 + croot * bw;
  
  acg.Sref = Sref;
  acg.bref = b;
  acg.cref = Sref / b;
  
  % LE points [x, y, z]
  % Root LE is now set to x=0. Tips are at x=0.0571 (swept/offset).
  
  p_tip_left  = [0.0571, -b/2,  0];
  p_root_left = [0.0, -bw/2, 0];
  p_root_right= [0.0,  bw/2, 0];
  p_tip_right = [0.0571,  b/2,  0];

  % --- Configuration ---
  nw = 60;           % COURSENESS: total number of vortex segments
  
  % MOMENT REFERENCE POINT: [x, y, z] 
  % Setting [0.0939 -0.0305 0.0] matches your balance point relative to Root LE (0,0,0).
  acg.pref = [0.0939 -0.0305 0.0]; 

  % --- Mesh Generation ---
  % Left Wing
  [acg, idx1] = add_wing(acg, p_tip_left, p_root_left, ctip, croot, 0, 0, floor(nw/3));
  
  % Mid Section (Fuselage/Rectangular)
  [acg, idx2] = add_wing(acg, p_root_left, p_root_right, croot, croot, 0, 0, floor(nw/3));
  
  % Right Wing
  [acg, idx3] = add_wing(acg, p_root_right, p_tip_right, croot, ctip, 0, 0, floor(nw/3));

  acg.vxwing = [idx1, idx2, idx3];

  % Generate airfoil data
  acg.foil = @finitewing_foils;
  acg.mwafgrid = build_finitewing(); % Main wing airfoils
  
  % No tail for this model
  acg.htpafgrid = [];
  acg.vtpafgrid = [];
  acg.vxhtp = [];
  acg.vxvtp = [];

  acg.ffsl = 0.0; % No fuselage drag model

end
