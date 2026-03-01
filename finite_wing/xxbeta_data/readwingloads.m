function [abq, bload, wload, coeff, span, wingarea, chord] = readwingloads(fname)
% Same as Draken using I676 balance 
%
% abq   : [alfa, beta, q] in columns (ncase x 3) 
% bload : body-axis loads (N m T C n l) 
% wload : wind-axis loads (D C L)
% coeff : wind-axis coefficients (CD CC CL Cl Cm Cn)
% span  : Model span [m]
% wingarea: Model reference area [m*m]
% chord : Model reference chord [m] 
% All moments are for the balance reference point.

  fid = fopen(fname,'r');  
  tem = fscanf(fid,'%f\n');
  fclose(fid);
  
  % not used
  span = tem(1);
  wingarea = tem(2);
  chord = tem(3);
  boffset = tem(4);
  
  % read balance calibration matrix
  % (not used in this version)
  k = 4;
  calib = zeros(40,6);
  for i=1:40
    for j=1:6
      k = k+1;
      calib(i,j) = tem(k);
    end
  end
  
  % Read loads (N m T C n l), alfa, voltages for tara (no wind, q=0)
  k = k+1;
  ntara=tem(k);
  for i=1:ntara
    for j=1:6
      k = k+1;
      tvolt(i,j) = tem(k);
    end
    k = k+1;
    talfa(i) = tem(k);
  end
  
  % read alfa beta q  
  % body_loads (N m T C n l)
  % wind_loads (D C L)
  % load_coefficients (CD CC CL Cl Cm Cn)
  % voltages (Ch1 Ch2 Ch3 Ch4 Ch5 Ch6)
  k = k+1;
  nrun = tem(k);
  for i=1:nrun
    k = k+1;
    alfa(i) = tem(k);
    k = k+1;
    beta(i) = tem(k);
    k = k+1;
    q(i) = tem(k);
    for j=1:6
      k = k+1;
      bload(i,j) = tem(k);
    end
    for j=1:3
      k = k+1;
      wload(i,j) = tem(k);
    end
    for j=1:6
      k = k+1;
      coeff(i,j) = tem(k);
    end
    for j=1:6
      k = k+1;
      volt(i,j) = tem(k);
    end
  end
  
  % assemble states
  abq = zeros(nrun, 3);
  abq(:,1) = alfa;
  abq(:,2) = beta;
  abq(:,3) = q;

