load detas

close all
h=axes
plot(v_236,de_236, "linewidth",3,'b-x', ...
     v_256,de_256, "linewidth",3,'r-x', ... 
     v_275,de_275, "linewidth",3,'g-x', ...
     v_294,de_294, "linewidth",3,'k-x')
axis([75 210 0 10])
hy=ylabel('Elevator \delta_e (deg)')
hx=xlabel('Airspeed (km/h)')
hl=legend('x_{cg} = 236 mm', ...
       'x_{cg} = 256 mm', ...
       'x_{cg} = 275 mm', ...
       'x_{cg} = 294 mm',"location",'northwest')
ifont=14
set(h,'FontSize',ifont)
set(hx,'FontSize',ifont)
set(hy,'FontSize',ifont)
set(hl,'FontSize',ifont)
