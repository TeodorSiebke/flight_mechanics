close all
if exist('t')
disp('t is set');
  else
    t=0:0.1:20;
end
e1 =  -0.3700 + 0.8973i;
e2 =   0.0034 + 0.0675i;
x=exp(e1*t)+0.001*exp(e2*t);
plot(t,real(x))
xlabel('t (s)');ylabel('Amplitude');
