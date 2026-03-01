for it=1:length(tvec)
x=[Y(it,:) interp1(fsm.de_set(:,1),fsm.de_set(:,2),tvec(it)) 0];
xdot=fplmod(tvec(it),x,fsm);
vtrue(it)=xdot(8);
nz(it)=xdot(13);
cl(it)=xdot(12);
alfdeg(it)=xdot(9)*180/pi;
end
close all
plot(Y(:,5),Y(:,6))
