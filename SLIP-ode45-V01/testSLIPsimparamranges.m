clear; close all; clc;

vbelt1 = 1.5/sqrt(9.81*0.88);
vbelt2 = 0.5/sqrt(9.81*0.88);
% vbelt1 = 1.5;
% vbelt2 = 0.5;

theta01 = -pi/6;
theta02 = -pi/4;
U01 = .5;
V01 = .6;
c1 = 0;
c2 = 0;
K1 = 1;
K2 = 1;

[ti,xi,yi,ui,vi,T1,T2,t2flight,d,discrep] = ...
    SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2);

plot(xi,yi,'b-');

% plotSLIPsim(ti,xi,yi,T1,T2,t2flight,vbelt1,vbelt2,d)