function [T,Z,te,ze,ie] = SLIPsim(U0,V0,theta0,K,c)
%% PROBLEM INFO, EQUATIONS OF MOTION
% This code solves SLIP from TOUCHDOWN to LIFTOFF
% Z = [x,y,u,v];
% x,u,y,v,theta
% L(t) = sqrt(X^2 + Y^2)
% dX = U;
% dU = K*(L0 - L)*sin(theta);
% sin(theta) = x/L;
% cos(theta) = y/L;
% dY = V;
% dV = K*(1-L)*cos(theta) - 1;
% theta is measured from the horizontal

%% PROBLEM PARAMETERS
% L0 = 1;
% g = 1;
% m = 1;
% These parameters normalize the problem and are hard-coded into the simulation. 

% ode45 event parameters
Tmax = 40;
f = 1.01; % this is for failure event1 (exceed leg length) in springmassTerminate.m

% Initial conditions (touchdown)

X0 = cos(theta0); %theta is positive CCW
Y0 = sin(theta0);
Z0 = [X0 Y0 U0 V0]; % define initial conditions

% this is the event function
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,f)); 

% this solves the ode system
[T,Z,te,ze,ie] = ode45(@(T,Z) springmassODE(T,Z,K,c),[0;Tmax],Z0,options);
