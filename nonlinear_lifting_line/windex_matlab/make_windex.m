function acg = make_windex()

% store reference properties for the Windex

  acg.Sref = 7.41;
  acg.bref = 12.1;
  acg.cref = 0.641;

% Windex wing geometry data
% Twist data uncertain, from SOR analysis sheets
% Root at 2.67 deg trailing edge down, then washout along span
% 0.5 deg at y=3.5 and 2 deg at tip
%
%       Stn 1     Stn 9     Stn 9     Stn 13    Stn 14

xwle = [0.0000    0.0500    0.0500    0.1300    0.1600];
ywle = [0.3000    3.5000    3.5000    5.5000    6.0000];
zwle = [0.0157    0.1832    0.1832    0.2878    0.3140];
wcrd = [0.7283    0.6416    0.6415    0.4000    0.3000];
twist= [2.6700    2.1700    2.1700    0.9700    0.6700]*pi/180;
       
% Horizontal tail plane (stabilizer)

xhtple=[2.7325    2.9215];
yhtple=[0.0000    1.0750];
zhtple=[0.9800    0.9800];
htpcrd=[0.5350    0.2650];

% Vertical tail plane (Fin)

xvtple=[1.6314    2.7325];
yvtple=[0.0000    0.0000];
zvtple=[0.0120    0.9800];
vtpcrd=[1.2018    0.5000];

% aerodynamic reference point; all moments computed about this pt.

  acg.pref = [0.0 0.0 0.0];

  % "mesh" : number of vortex segments for main wing and tail
  nw = 28; 
  nt = nw/2;
  nvtp = 10;
  acg.nw=nw;
  
% left wing (Note, discontinuity in data at midwing station)
	     
  [acg, idx1] = add_wing(acg, [xwle(5), -ywle(5), zwle(5)], ...
		      [xwle(4), -ywle(4), zwle(4)], ...
                 wcrd(5), wcrd(4), twist(5), twist(4), 4);
 
  [acg, idx2] = add_wing(acg, [xwle(4), -ywle(4), zwle(4)], ...
		      [xwle(3), -ywle(3), zwle(3)], ...
                 wcrd(4), wcrd(3), twist(4), twist(3), 4);
 
  [acg, idx3] = add_wing(acg, [xwle(2), -ywle(2), zwle(2)], ...
		      [xwle(1), -ywle(1), zwle(1)], ...
                 wcrd(2), wcrd(1), twist(2), twist(1), 4);
 
% Fuselage section

  [acg, idx4] = add_wing(acg, [xwle(1), -ywle(1), zwle(1)], ...
                      [xwle(1),  ywle(1), zwle(1)], ...
                 wcrd(1), wcrd(1), twist(1), twist(1), 4);

% right wing

  [acg, idx5] = add_wing(acg, [xwle(1),  ywle(1), zwle(1)], ...
		      [xwle(2), ywle(2), zwle(2)], ...
                 wcrd(1), wcrd(2), twist(1), twist(2), 4);

  [acg, idx6] = add_wing(acg, [xwle(3),  ywle(3), zwle(3)], ...
		      [xwle(4), ywle(4), zwle(4)], ...
                 wcrd(3), wcrd(4), twist(3), twist(4), 4);

  [acg, idx7] = add_wing(acg, [xwle(4),  ywle(4), zwle(4)], ...
		      [xwle(5), ywle(5), zwle(5)], ...
                 wcrd(4), wcrd(5), twist(4), twist(5), 4);

  acg.vxwing = [idx1, idx2, idx3, idx4, idx5, idx6, idx7];
 
  % Flaps: y=0.3 to y=3.5 (idx3, idx5) + Fuselage (idx4)
  % Ailerons: y=3.5 to y=6.0 (idx2+idx1, idx6+idx7)
  
  acg.vxflap = [idx3, idx4, idx5]; 
  acg.vxail_l = [idx1, idx2];      % Left Aileron (Outer sections)
  acg.vxail_r = [idx6, idx7];      % Right Aileron (Outer sections)


% HTP - in neutral position 

  thtp=0.0*pi/180;
  [acg, idxl] = add_wing(acg, [xhtple(2), -yhtple(2) zhtple(2)], ...
                 [xhtple(1), yhtple(1), zhtple(1)], ...
                 htpcrd(2), htpcrd(1), thtp, thtp, nt/2);
  [acg, idxr] = add_wing(acg, [xhtple(1), yhtple(1) zhtple(1)], ...
                 [xhtple(2), yhtple(2), zhtple(2)], ...
                 htpcrd(1), htpcrd(2), thtp, thtp, nt/2);
  acg.vxhtp = [idxl, idxr];

% Model vertical tail, the fin

  [acg, idx] = add_wing(acg,[xvtple(1),yvtple(1),zvtple(1)], ...
                     [xvtple(2),yvtple(2),zvtple(2)], ...
                     vtpcrd(1),vtpcrd(2),0.0,0.0,nvtp);
  acg.vxvtp = idx;

% Generate airfoil data

  acg.foil = @windex_foils;
  acg.mwafgrid=build_wxwing();   % Main wing airfoils
  acg.htpafgrid=build_wxtail();  % HTP airfoils
  acg.vtpafgrid=build_wxtail();  % VTP airfoils
  
% fuselage drag in the form S*CD
  % Isolated fuselage drag from wind tunnel u=35m/s: CD_exp - CD_sim = 0.00748 - 0.00533 = 0.00215
  acg.ffsl = 0.00215 * acg.Sref;

end 

