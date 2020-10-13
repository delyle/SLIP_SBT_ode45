K = 10;
c = 1;
theta0 = pi/2;
vbelt = 0;
U0 = -0.4+vbelt;
V0 = -0.1;


[T,Z,te,ze,ie] = SLIPsim(U0,V0,theta0,K,c);

close all
plotSLIPsimOnePhase(T,Z,c,K,vbelt)