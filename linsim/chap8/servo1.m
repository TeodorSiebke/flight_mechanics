% First order servo model

tau=1
A=[-1/tau];B=[1/tau];
C=[1];D=[0];

close all

sys=ss(A,B,C,D);
step(sys,10)
