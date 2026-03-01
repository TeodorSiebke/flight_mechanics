function acg = make_windexfusewing()

% Span for this analysis is actual - fuselage width: 600-100=500 mm
% store reference properties for the Windex, scale 1:5

  acg.Sref = 7.41/25.0;
  acg.bref = 12.1/5.0;
  acg.cref = 0.641/5.0;

% Windex wing geometry data
% Root at 2.67 deg trailing edge down
%       root      tip

xwle = [0.0000    0.0];
ywle = [0.0000    0.3]; % Semispan 250mm
zwle = [0.0       0.0];
wcrd = [0.1385    0.1385];
twist= [2.6700    2.6700]*pi/180;
       
% aerodynamic reference point; all moments computed about this pt.
% This balance reference point

  acg.pref = [0.0354 0.0 0.0];

  % "mesh" : number of vortex segments for main wing and tail
  nw = 10; 
  acg.nw=nw;
  
% left wing (Note, discontinuity in data at midwing station)
	     
% Fuselage section

  [acg, idx1] = add_wing(acg, [xwle(2), -ywle(2), zwle(2)], ...
                      [xwle(2),  ywle(2), zwle(2)], ...
                 wcrd(1), wcrd(1), twist(1), twist(1), nw);

  acg.vxwing = [idx1];

% Generate airfoil data

  acg.foil = @windexfusewing_foils;
  acg.mwafgrid=build_wxwing();   % Main wing airfoils
  
% fuselage drag in the form S*CD

  acg.ffsl = 0.0;

end 

