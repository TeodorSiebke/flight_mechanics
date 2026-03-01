function [cd,cl,cm]=lon_aero(alpha,deltae,qbar,fsm)

% alfa and deltae in radians
% Simple model with linear approx in deltae for CL and Cm
% Assume drag independent of deltae (linear approximation is useless)

  if alpha > fsm.valpha(end), error('alpha too large'), end
  if alpha < fsm.valpha(1), error('alpha too small'), end

  
vcl=interp1(fsm.valpha,fsm.cldat,alpha);
vcm=interp1(fsm.valpha,fsm.cmdat,alpha);
vcd=interp1(fsm.valpha,fsm.cddat,alpha);

clde=(vcl(3)-vcl(1))/(fsm.vde(3)-fsm.vde(1));
cl=vcl(2)+clde*deltae+fsm.CLq*qbar;

cmde=(vcm(3)-vcm(1))/(fsm.vde(3)-fsm.vde(1));
cm=vcm(2)+cmde*deltae+fsm.Cmq*qbar;

cd=vcd(2);

