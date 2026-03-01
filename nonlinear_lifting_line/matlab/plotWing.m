function [] = plotWing(acg,nw,nw_t)
%Plots the wing geometry
%nw is amount of vortex elements in main wing
%nw_t is amount of vortex elements in tailplane

%Main Wing
qC = [acg.pa(1,1:2); 0.5*(acg.pa(1:nw,1:2)+acg.pb(1:nw,1:2)); acg.pb(nw,1:2)];  %Quarter Chord
LE = [qC(:,1) - 0.25*[0.2134; acg.chord(1:nw); 0.2134],qC(:,2)] ;    %Leading Edge
TE = [qC(:,1) + 0.75*[0.2134; acg.chord(1:nw); 0.2134],qC(:,2)] ;    %Trailing Edge
FH = [qC(:,1) + 0.60*[0.2134; acg.chord(1:nw); 0.2134],qC(:,2)] ;    %Flap Hinge

%Plot Main Wing
figure
plot(FH(:,1),FH(:,2),'k--','LineWidth',1.5)   %FLHI

hold on
plot(LE(:,1),LE(:,2),'b','LineWidth',1.5)
plot(TE(:,1),TE(:,2),'b','LineWidth',1.5)
plot([LE(1,1) TE(1,1)],[LE(1,2) TE(1,2)],'b','LineWidth',1.5)  
plot([LE(end,1) TE(end,1)],[LE(end,2) TE(end,2)],'b','LineWidth',1.5)
view([90 90])

%Tail Plane
qC_t = [acg.pa(nw+1,1:2); 0.5*(acg.pa(nw+1:end,1:2)+acg.pb(nw+1:end,1:2)); acg.pb(end,1:2)];
LE_t = [qC_t(:,1) - 0.25*[0; acg.chord(nw+1:end); 0],qC_t(:,2)] ;
TE_t = [qC_t(:,1) + 0.75*[0; acg.chord(nw+1:end); 0],qC_t(:,2)] ;
FH_t = [qC_t(:,1) + 0.25*[0; acg.chord(nw+1:end); 0],qC_t(:,2)] ;    %Flap Hinge

%Plot HTP
%plot(FH_t(:,1),FH_t(:,2),'k--')   %25% Chord
hold on
plot(LE_t(:,1),LE_t(:,2),'b','LineWidth',1.5)
plot(TE_t(:,1),TE_t(:,2),'b','LineWidth',1.5)
plot([LE_t(1,1) TE_t(1,1)],[LE_t(1,2) TE_t(1,2)],'b','LineWidth',1.5)  
plot([LE_t(end,1) TE_t(end,1)],[LE_t(end,2) TE_t(end,2)],'b','LineWidth',1.5)
view([90 90])

%plot([0.7 4.2],[-0.25 -0.25],'b');
%plot([0.7 4.2],[0.25 0.25],'b');
axis equal
axis([-0.5 5 -10 10])
set(gca,'FontSize',20)
legend('Flap Hinge at 85 % of chord',[330 110 1 1]);
ylabel('y [m]')
xlabel('x [m]')
