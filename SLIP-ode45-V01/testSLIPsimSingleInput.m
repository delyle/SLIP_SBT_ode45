% This is script is meant for testing purposes (i.e. to run SLIPsim.m
% and determine the permissible ranges of inputs that make sense for the
% non-dimensionalization of this problem.
clear; close all; clc;

%% Non-dimensionalization of problem

% [M] = m
% [L] = L0
% [T] = sqrt(L0/g)

% X' = x/L0
% U' = u/sqrt(g*L0)
% K' = k/(m*g/L0)
% C' = c/(m*sqrt(g/L0))

m = 80; %in kg
L0 = 1; %in m
g = 9.81; %in m/s^2

%fully animate the solution, distinguish phase 1/2

%% Simulation
% all input values are non-dimensional by default

%Read in values
vbelt1 = 1.5/sqrt(g*L0);
vbelt2 = 0.5/sqrt(g*L0);
maxBeltSpd = max(vbelt1,vbelt2);
% theta is positive CCW, where zero is the vertical
theta01 = pi/3;
theta02 = pi/3;
U01 = 0.5*maxBeltSpd;
V01 = 0.5*maxBeltSpd;

K = 10;
K1 = K;
K2 = K;

% zeta range of interest
% zeta = linspace(0,2,5);
zeta = 0;

c = zeta*(2*sqrt(K)); % calculate C' based on zeta for input into the problem
c1 = c;
c2 = c;

[ti,xi,yi,ui,vi,T1,T2,t2flight,d,discrep] = ...
                    SLIPsim(vbelt1,vbelt2,theta01(i),theta02(j),U01(k),V01(l),c1,c2,K1,K2);
                if isnan(discrep)
                    disp('solution invalid');
                end

                if ~isnan(discrep)
                    SLIPanim(ti,xi,yi,ui,vi,T1,T2,t2flight,d,vbelt1,vbelt2);
                end

