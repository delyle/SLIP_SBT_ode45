clear; close all; clc;

m = 80; %in kg
L0 = 1; %in m
g = 9.81; %in m/s^2

%fully animate the solution, distinguish phase 1/2

vbelt1 = 1.5/sqrt(g*L0);
vbelt2 = 0.5/sqrt(g*L0);
% vbelt1 = 1.5;
% vbelt2 = 0.5;

theta01 = -pi/4;
theta02 = -pi/4;
U01 = 1;
V01 = .5;
c1 = 0;
c2 = 0;

K = 10;
K1 = K;
K2 = K;

[ti,xi,yi,ui,vi,T1,T2,t2flight,d,discrep] = ...
    SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2);

plot(xi,yi,'b-');

% plotSLIPsim(ti,xi,yi,T1,T2,t2flight,vbelt1,vbelt2,d)