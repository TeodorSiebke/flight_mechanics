function [xdot]=fplmod(tval,xvar,fsm)

% This version with both state and control variables in xvar
% Returned functions are both state and output equations
% Use a wrapper for most types of analysis

% Some constants in use

zero=0.0;half=0.5;one=1.0;two=2.0;
r2d=57.29577951308232;kappa=1.4;
expg1=(kappa-one)/kappa;
expg2=kappa/(kappa-one);
aspeed0=340.269;press0=101325.0;
tusen=1000.0;

% Aircraft data

      cbar = fsm.cref;         % Reference chord (m)
      span = fsm.bref;         % Span (m)
      Sref = fsm.Sref;         % Reference area (m*m)
      gacc = 9.81;             % Acceleration of gravity (m/(s*s))
      tepsr= 0.0*pi/180;       % Thrust installation angle in radians

%     Define all the variables

      u        = xvar(1);    % m/s
      w        = xvar(2);    % m/s
      qrate    = xvar(3);    % rad/s
      theta    = xvar(4);    % rad
      xdiste   = xvar(5);    % m
      altm     = xvar(6);    % m
      fuelmass = xvar(7);    % kg

%     Control variables 

      deltae   = xvar(8);    % rad
      deltap   = xvar(9);    % 0 (flight idle) to 1 (full thrust) 

%     Aircraft mass

      amass=fsm.mass;

%     Reference

      xref=fsm.pref(1);
      yref=fsm.pref(2);
      zref=fsm.pref(3);

      if isfield(fsm,'cog')
        % Using the shift logic: x_aft (Aero) to x_fwd (FM) and z_up (Aero) to z_down (FM)
        xcg=fsm.cog(1); ycg=fsm.cog(2); zcg=fsm.cog(3);
      else
        xcg=xref;ycg=yref;zcg=zref;
      end

%     Transformation matrix body to earth

      sintet=sin(theta);
      costet=cos(theta);

%     Compute wind axis state

      airspd = sqrt( u*u + w*w );
      alpha   = atan(w/u);

%     More extra variables

      sinalf=sin(alpha);
      cosalf=cos(alpha);

% Pitch mass moment of inertia

      aIy = fsm.massmom(2,2);       % Pitching moment of inertia

%     Define some temporary variables

      altkm=altm/1000;       % Altitude in km

%     Get atmospheric data
      alt_safe = max(0, min(altm, 50000));
      [rho,aspeed,temp,press] = isaatm(alt_safe);

%     Compute Mach number

      rmach=airspd/aspeed;
      rmach2=rmach*rmach;
      qdyn=rho*airspd*airspd/two;

%     Compute calibrated (indicated) air speed
%     First, compute the total pressure that the aircraft should sense

      ptot=press*(one+((kappa-one)/two)*rmach2)^expg2;
      qc=ptot-press;
      rmtem=sqrt( (two/(kappa-one))*( (qc/press0 + one)^expg1 - one ) );
 
%     Indicated (calibrated) airspeed and derivatives

      vcal=rmtem*aspeed0;

%     Nondimensional pitch rate

      qbar=qrate*cbar/(2*airspd);

%     Engine data (Simple scaling of full thrust)

      thrust=fsm.thrustmax*deltap;
      fuelb=0.0;  % Electric

%     Longitudinal aerodynamics interpolation
      [CD,CL,Cm]=lon_aero(alpha,deltae,qbar,fsm,airspd);

      % Add Flap effects (assuming dCL_df, dCD_df, dCm_df are per rad)
      % Use xvar(11) as a temporary place for flaps if needed, or pass via fsm
      if length(xvar) >= 11
          deltaf = xvar(11);
      else
          deltaf = 0;
      end
      
      CL = CL + fsm.dCL_df * deltaf;
      CD = CD + fsm.dCD_df * deltaf;
      Cm = Cm + fsm.dCm_df * deltaf;

      rlift = qdyn*Sref*CL;

      drag = qdyn * Sref *CD;

      Cmadot = fsm.Cmadot;

%     Stevens and Lewis model as on p. 123-124.
%     Forces in body axis system

      xforce = thrust - drag*cosalf  + rlift*sinalf;
      zforce =        - rlift*cosalf - drag*sinalf;

%     Force equation(s) body axis

      udot = - qrate*w + xforce/amass - gacc*sintet;
      wdot = qrate*u + zforce/amass + gacc*costet;

%     Time derivative of airspeed

      aspdot = (u*udot + w*wdot)/airspd;

%     Time derivative of angle-of-attack

      alfdot = (cosalf*cosalf)*(wdot*u-w*udot)/(u*u);

%     Moment equation(s)

pitmom = (qdyn*Sref*cbar*(Cm + Cmadot*alfdot*cbar/(2*airspd))+thrust*fsm.z_eng);

%     Implement adjust for c.g. position (using Aero coordinates for xcg, xref, zcg, zref)
%     Since M = z_pos*X_force - x_pos*Z_force in FM frame:
shift_x = xcg - xref; % Positive if CG is aft of ARP
shift_z = zcg - zref; % Positive if CG is above ARP

pitmom = pitmom + shift_z * xforce - shift_x * zforce;

%     Vehicle symmetry is assumed

      qdot = ( pitmom )/aIy;

%     Where we are equations

      xedot  =  costet*u + sintet*w;
      altdot = - ( -sintet*u + costet*w);

      tetdot = qrate ;

%     State equations

      xdot(1:13,1)=0;

      xdot(1) = udot;
      xdot(2) = wdot;
      xdot(3) = qdot;
      xdot(4) = tetdot;
      xdot(5) = xedot;
      xdot(6) = altdot;
      xdot(7) = -fuelb;

%     Some useful additional information

      % Magnitude of all g force vectors in all directions
      nz = sqrt(xforce^2 + zforce^2) / (amass * gacc);

%     Extra output variables, add more if you like

      xdot(8)=airspd;
      xdot(9)=alpha;
      xdot(10)=vcal;
      xdot(11)=rmach;
      xdot(12)=CL;
      xdot(13)=nz;

