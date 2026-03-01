function [xdot]=simsub(t,x,fsm)

% Wrapper function for ODE integration

xs=zeros(9,1); % Initialize all variables to zero

xs(1)=x(1);     % u
xs(2)=x(2);     % w
xs(3)=x(3);     % q
xs(4)=x(4);     % theta
xs(5)=x(5);     % xe
xs(6)=x(6);     % alt
xs(7)=x(7);     % fuelmass

xs(8)=interp1(fsm.de_set(:,1),fsm.de_set(:,2),t);    % deltae
xs(9)=fsm.dp_set;    % deltap

% Evaluate all state functions

fs=fplmod(0,xs,fsm);

% Define the differential equations

xdot(1,1)=fs(1);
xdot(2,1)=fs(2);
xdot(3,1)=fs(3);
xdot(4,1)=fs(4);
xdot(5,1)=fs(5);
xdot(6,1)=fs(6);
xdot(7,1)=fs(7);
