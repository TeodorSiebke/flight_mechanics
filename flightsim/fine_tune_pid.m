% fine_tune_pid.m
addpath('flightsim');
v_target=25; h_target=0; fsm=make_fsim(); 
g0=[0.0045, 0.0002, 0.015, 0.4]; 
options=optimset('Display','iter','TolX',1e-5); 
[best_g, fval] = fminsearch(@(g) calculate_point_cost(g, v_target, h_target, fsm), g0, options); 
disp('Best Gains:'); disp(best_g);

function c = calculate_point_cost(g, V, H, fsm)
    Kp=g(1); Ki=g(2); Kd=g(3); Kq=g(4); 
    if any(g<=0); c=1000; return; end
    x_tr=[V 0.5 0 0.05 0 H 0 0.06 0.1]'; 
    ivar=[1 2 4 8 9]; ifun=[1 2 3 6 8]; xt=x_tr(ivar); 
    for iter=1:30; 
        x_tr(ivar)=xt; xd=fplmod(0,x_tr,fsm); 
        res=xd(ifun); res(5)=res(5)-V; 
        if norm(res)<1e-8, break; end; 
        xs=1e-7; J=zeros(5,5); 
        for k=1:5; 
            xh=x_tr; xh(ivar(k))=xh(ivar(k))+xs; xmh=x_tr; xmh(ivar(k))=xmh(ivar(k))-xs; 
            J(:,k)=(fplmod(0,xh,fsm)-fplmod(0,xmh,fsm))/(2*xs); 
        end; 
        xt = xt - 0.7 * (J(ifun,:)\res); 
    end; 
    x_tr(ivar)=xt; 
    s_idx=[1 2 3 4 6]; n=6; A=zeros(n,n); xsl=1e-5; 
    cls=@(xv) [fplmod(0,[xv(1);xv(2);xv(3);xv(4);0;xv(5);0;x_tr(8)-(Kp*(H-xv(5))+Ki*xv(6)-Kd*(-(-sin(xv(4))*xv(1)+cos(xv(4))*xv(2))))+Kq*xv(3);x_tr(9)],fsm); H-xv(5)]; 
    for k=1:n; 
        xh=[x_tr(s_idx);0]; xh(k)=xh(k)+xsl; xmh=[x_tr(s_idx);0]; xmh(k)=xmh(k)-xsl; 
        xdh=cls(xh); xdh=[xdh(1:4);xdh(6);xdh(14)]; 
        xdm=cls(xmh); xdm=[xdm(1:4);xdm(6);xdm(14)]; 
        A(:,k)=(xdh-xdm)/(2*xsl); 
    end; 
    eg=eig(A); 
    c = max(real(eg)); 
end
