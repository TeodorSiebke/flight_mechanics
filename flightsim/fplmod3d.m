function [xdot, extra]=fplmod3d(tval,xstate,fsm)

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

%     Aircraft mass

      amass=fsm.mass;

%     Define all the variables

      u      = xstate(1);    % m/s
      v      = xstate(2);    % m/s
      w      = xstate(3);    % m/s
      prate  = xstate(4);    % rad/s
      qrate  = xstate(5);    % rad/s
      rrate  = xstate(6);    % rad/s
      phi    = xstate(7);    % rad
      theta  = xstate(8);    % rad
      psi    = xstate(9);    % rad
      xdiste = xstate(10);   % m
      ydiste = xstate(11);   % m
      altm   = xstate(12);   % m

      deltae = xstate(13);   % rad
      deltaa = xstate(14);   % rad
      deltar = xstate(15);   % rad
      deltap = xstate(16);   % rad

      altint = 0;
      velint = 0;
      drdot  = 0;
      psiint = 0;

% altint   is additional state for integrator term (altitude)
% velint   is additional state for integrator term (velocity)
% psiint   is additional state for integrator term (heading)

%     Extra intermediate variables

sinphi=sin(phi);
cosphi=cos(phi);

sintet=sin(theta);
costet=cos(theta);
tantet=tan(theta);
sectet=one/costet;

sinpsi=sin(psi);
cospsi=cos(psi);

%     Transformation matrix body to earth

tvbmx(1,1)=costet*cospsi;
tvbmx(1,2)=sinphi*sintet*cospsi-cosphi*sinpsi;
tvbmx(1,3)=cosphi*sintet*cospsi+sinphi*sinpsi;
tvbmx(2,1)=costet*sinpsi;
tvbmx(2,2)=sinphi*sintet*sinpsi+cosphi*cospsi;
tvbmx(2,3)=cosphi*sintet*sinpsi-sinphi*cospsi;
tvbmx(3,1)=-sintet;
tvbmx(3,2)=sinphi*costet;
tvbmx(3,3)=cosphi*costet;

%     Wind blowing

xwdot=zero;
ywdot=zero;
zwdot=zero;
xwind=zero;
ywind=zero;
zwind=zero;

%     Compute wind axis state

vspeed = sqrt( u*u + v*v + w*w );

uair = u - tvbmx(1,1)*xwind - tvbmx(2,1)*ywind -tvbmx(3,1)*zwind;
vair = v - tvbmx(1,2)*xwind - tvbmx(2,2)*ywind -tvbmx(3,2)*zwind;
wair = w - tvbmx(1,3)*xwind - tvbmx(2,3)*ywind -tvbmx(3,3)*zwind;

airspd = sqrt( uair*uair + vair*vair + wair*wair );
alpha  = atan(wair/uair);
beta   = asin(vair/airspd);

%     More extra variables

sinalf=sin(alpha);
cosalf=cos(alpha);

sinbet=sin(beta);
cosbet=cos(beta);
tanbet=tan(beta);
secbet=one/cosbet;

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

      z_eng=fsm.z_eng;                 % Engine c.g. offset

      aIx    = fsm.massmom(1,1);
      aIy    = fsm.massmom(2,2);       % Pitching moment of inertia
      aIz    = fsm.massmom(3,3);
      aIzx   = fsm.massmom(3,1);

%     Define some temporary variables

      altkm=altm/1000;       % Altitude in km

%     Get atmospheric data

      alt_safe = max(0, min(altm, 50000));
      [rho,aspeed,temp,press] = isaatm(alt_safe);

%     Correct for wind here

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

%     Aerodynamics

      pbar=prate*span/(2*airspd);
      qbar=qrate*cbar/(2*airspd);
      rbar=rrate*span/(2*airspd);

%     Engine data

      thrust=fsm.thrustmax*deltap;
      fuelb=0.0;  % Electric

%     Interpolate functions for aerodata

      [CD,CL,Cm]=lon_aero(alpha,deltae,qbar,fsm);

      % Add Flap effects (assuming dCL_df, dCD_df, dCm_df are per rad)
      % We need a flap deflection. Let's look for it in xstate(17) or fixed.
      % The prompt mentions "flap and aileron control". 
      % Let's assume xstate(17) is flap.
      if length(xstate) >= 17
          deltaf = xstate(17);
      else
          deltaf = 0;
      end
      
      CL = CL + fsm.dCL_df * deltaf;
      CD = CD + fsm.dCD_df * deltaf;
      Cm = Cm + fsm.dCm_df * deltaf;

      rlift = qdyn*Sref*CL;

%     Drag

      drag = qdyn*Sref*CD;

%     Moment around aerodynamic reference point

      Cmadot = fsm.Cmadot;

%     Lateral ADD CONTROL

      yforce = fsm.Cybeta*beta + fsm.Cyp*pbar + fsm.Cyr*rbar + ...
               fsm.Cyda*deltaa + fsm.Cydr*deltar;

      yforce = qdyn*Sref*yforce;

      cl = fsm.Clbeta*beta + fsm.Clp*pbar + fsm.Clr*rbar + ...
           fsm.Clda*deltaa + fsm.Cldr*deltar;

      cn = fsm.Cnbeta*beta + fsm.Cnp*pbar + fsm.Cnr*rbar + ...
           fsm.Cnda*deltaa + fsm.Cndr*deltar;

%     Stevens and Lewis model as on p. 123-124.

%     Forces in body axis system

      xforce = thrust - drag*cosalf  + rlift*sinalf;
      zforce =        - rlift*cosalf - drag*sinalf;

%     Force equation(s) body axis

      udot = rrate*v - qrate*w + xforce/amass - gacc*sintet;
      vdot = prate*w - rrate*u + yforce/amass + gacc*costet*sinphi;
      wdot = qrate*u - prate*v + zforce/amass + gacc*costet*cosphi;

%     Time derivative of ground speed

      veldot = (u*udot + v*vdot + w*wdot)/vspeed;

%     Time derivative of airspeed

      uadot = udot - rrate*ywind + qrate*zwind ...
              - tvbmx(1,1)*xwdot - tvbmx(2,1)*ywdot -tvbmx(3,1)*zwdot;

      vadot = vdot + rrate*xwind - prate*zwind ...
              - tvbmx(1,2)*xwdot - tvbmx(2,2)*ywdot -tvbmx(3,2)*zwdot;

      wadot = wdot - qrate*xwind + prate*ywind ...
              - tvbmx(1,3)*xwdot - tvbmx(2,3)*ywdot -tvbmx(3,3)*zwdot;

      aspdot = (uair*uadot + vair*vadot + wair*wadot)/airspd;

%     Time derivative of angle-of-attack

      alfdot = (one+alpha*alpha)*(wadot*uair-wair*uadot)/(uair*uair);

%     Moment equation(s)

pitmom = (qdyn*Sref*cbar*(Cm + Cmadot*alfdot*cbar/(2*airspd))+thrust*z_eng);

%     Adjust for c.g. position (using Aero coordinates for xcg, xref, zcg, zref)
shift_x = xcg - xref; 
shift_z = zcg - zref; 

pitmom = pitmom + shift_z * xforce - shift_x * zforce;

rolmom = qdyn*Sref*( cl ) * span;
yawmom = qdyn*Sref*( cn ) * span;

%     Vehicle symmetry is assumed

      qdot = ( pitmom + aIzx*(rrate*rrate-prate*prate) ...
                      + (aIz-aIx)*rrate*prate )/aIy;

%     Inverse of remaining inertia matrix is applied to obtain rotation rates

      detain=aIx*aIz - aIzx*aIzx;

      rhs(1) = rolmom + aIzx*prate*qrate + (aIy-aIz)*qrate*rrate ;
      rhs(2) = yawmom - aIzx*qrate*rrate + (aIx-aIy)*prate*qrate;

      pdot = (aIz*rhs(1)+aIzx*rhs(2))/detain;

      rdot = (aIzx*rhs(1)+aIx*rhs(2))/detain;

%     Where we are equations

      xedot  =  tvbmx(1,1)*u + tvbmx(1,2)*v + tvbmx(1,3)*w;
      yedot  =  tvbmx(2,1)*u + tvbmx(2,2)*v + tvbmx(2,3)*w;
      altdot = -tvbmx(3,1)*u - tvbmx(3,2)*v - tvbmx(3,3)*w;

      phidot = prate + (qrate*sinphi+rrate*cosphi) * tantet;
      tetdot = qrate - rrate*sinphi;
      psidot = (qrate*sinphi + rrate*cosphi)*sectet;

%     Wind system angles

      thetaw=asin(altdot/vspeed);
      costew=cos(thetaw);

%     acos() is dangerous, could be undefined if argument is 1+epsmch

      tem=min(one,yedot/(vspeed*costew));
      psiw=asin(tem);

%     State equations

      xdot = zeros(13,1);

      xdot(1) = udot; % veldot
      xdot(2) = vdot; % alfdot
      xdot(3) = wdot; % betdot
      xdot(4) = pdot;
      xdot(5) = qdot;
      xdot(6) = rdot;
      xdot(7) = phidot;
      xdot(8) = tetdot;
      xdot(9) = psidot;
      xdot(10) = xedot;
      xdot(11) = yedot;
      xdot(12) = altdot;
%      xdot(13) = -fuelb;

%     Some useful additional information

      extra.alpha = alpha * 180/pi;
      extra.nz = qdyn*Sref*CL/(amass*gacc);
      extra.CL = CL;

