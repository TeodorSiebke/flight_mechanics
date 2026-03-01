function [cd,cl,cm]=lon_aero(alpha,deltae,qbar,fsm,V)

% V in m/s (optional, for speed-dependent aero)
if nargin < 5
    V = 30;
end

% Clamping alpha
if alpha > fsm.valpha(end), alpha=fsm.valpha(end); end
if alpha < fsm.valpha(1), alpha = fsm.valpha(1); end

% Determine weighting for speed interpolation
v_anchors = fsm.v_anchors;
if V <= v_anchors(1)
    w_v = [1 0 0];
elseif V >= v_anchors(3)
    w_v = [0 0 1];
elseif V < v_anchors(2)
    frac = (V - v_anchors(1)) / (v_anchors(2) - v_anchors(1));
    w_v = [(1-frac) frac 0];
else
    frac = (V - v_anchors(2)) / (v_anchors(3) - v_anchors(2));
    w_v = [0 (1-frac) frac];
end

% Interpolate tables
cldat = w_v(1)*fsm.cldat_30 + w_v(2)*fsm.cldat_50 + w_v(3)*fsm.cldat_70;
cmdat = w_v(1)*fsm.cmdat_30 + w_v(2)*fsm.cmdat_50 + w_v(3)*fsm.cmdat_70;
cddat = w_v(1)*fsm.cddat_30 + w_v(2)*fsm.cddat_50 + w_v(3)*fsm.cddat_70;

% Interpolate for alpha
vcl = interp1(fsm.valpha, cldat, alpha);
vcm = interp1(fsm.valpha, cmdat, alpha);
vcd = interp1(fsm.valpha, cddat, alpha);

% Elevator derivatives (assume linear between anchor tables)
clde = (vcl(3)-vcl(1))/(fsm.vde(3)-fsm.vde(1));
cl = vcl(2) + clde*deltae + fsm.CLq*qbar;

cmde = (vcm(3)-vcm(1))/(fsm.vde(3)-fsm.vde(1));
cm = vcm(2) + cmde*deltae + fsm.Cmq*qbar;

cd = vcd(2);

