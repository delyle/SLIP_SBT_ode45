function [ti,xi,yi,ui,vi,L0,g,T1,T2,t2flight,d,Ef,E0,xprev,yprev] = SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,h01,c,K1,K2)
% This is the function run by SLIP_GUI to solve SLIP on a treadmill
% Function outputs are for use by the plotting function, plotSLIPsim.m
% For more details on how this function is used, see SLIP_GUI.mlapp

if nargin < 7
    K1 = 1;
    K2 = 1;
elseif nargin < 8
    K2 = K1;
end
    
%% INFO
% This code solves SLIP on a treadmill
% Z = [x,y,u,v];
% x,u,y,v,theta
% L(t) = sqrt(X^2 + Y^2)
% dX = U;
% dU = K*(L0 - L)*sin(theta);
% sin(theta) = x/L;
% cos(theta) = y/L;
% dY = V;
% dV = K*(1-L)*cos(theta) - 1;
% close all; clc; clear;
saveanimation = true;
filename = "SLIP.gif";

%% PROBLEM DATA
g = 1;
% vbelt1 = 1.5; % Belt speed (leftward)
% vbelt2 = 0.5;
% E0 = 20; % initial energy; constrains V0
Tmax = 40;
L0 = 10;
f = 1.01; %this is for the event1 (mass falling through floor)

%% PHASE 1
% theta01 = -pi/5; % angle of attack (radians, measured to vertical)
X01 = L0*sin(theta01);
Y01 = L0*cos(theta01);
% U01 = 2; % landing horizontal speed (rightward)
V01 = sqrt(2*(h01 - Y01));
%Check if Vland is valid
Z01 = [X01 Y01 U01 -V01];
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt1,f)); % this is the event function
[T1,Z1,te1,ze1,ie1] = ode45(@(T,Z) springmassODE(T,Z,K1,L0,vbelt1,c),[0;Tmax],Z01,options); % this solves the ode
X1 = Z1(:,1); Y1 = Z1(:,2); U1 = Z1(:,3); V1 = Z1(:,4);

%% PHASE 2
% theta02 = pi/2.6;
X02 = L0*sin(theta02);
Y02 = L0*cos(theta02);
U02 = U1(end); % U is conserved
Yapex1 = Y1(end) + V1(end)^2/(2*g);
V02 = sqrt(2*g*abs(Y02 - Yapex1)); %make negative since we know it is downwards
Z02 = [X02 Y02 U02 -V02];
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt2,f)); % this is the event function
[T2,Z2,te2,ze2,ie2] = ode45(@(T,Z) springmassODE(T,Z,K2,L0,vbelt2,c),[0;Tmax],Z02,options); % this solves the ode
X2 = Z2(:,1); Y2 = Z2(:,2); U2 = Z2(:,3); V2 = Z2(:,4);

%% JOIN INTO CONTINUOUS VECTOR
Tq1 = linspace(0,T1(end),101);
Tq2 = linspace(0,T2(end),101);
Zq1 = interp1(T1,Z1,Tq1,'spline');
Zq2 = interp1(T2,Z2,Tq2,'spline');


%Flight phase prior to Phase 1
Tprev = abs(V01/g);
tprev = linspace(0,Tprev);
xprev = X01 - U01*tprev(end) + U01.*tprev;
yprev = h01 - 1/2*g*tprev.^2;

%Redefine time and states with query values
T1 = Tq1'; T2 = Tq2';
X1 = Zq1(:,1); Y1 = Zq1(:,2); U1 = Zq1(:,3); V1 = Zq1(:,4);
X2 = Zq2(:,1); Y2 = Zq2(:,2); U2 = Zq2(:,3); V2 = Zq2(:,4);

%Flight Phase 1
T1apex = V1(end)/g; %time to reach apex
T1drop = abs(2*(Y02 - Yapex1)/abs(V02)); %time to drop from apex
T1flight = T1apex + T1drop;

d = X1(end) - X02 + U1(end)*T1flight; %"step length" to reset the origin

t1flight = linspace(0,T1flight,51)' + T1(end);
x1flight = X1(end) + U1(end)*(t1flight - T1(end));
y1flight = Y1(end) - 1/2*g*(t1flight - T1(end)).^2 + V1(end)*(t1flight-T1(end));
u1flight = U1(end)*ones(size(t1flight));
v1flight = V1(end) - g*(t1flight - T1(end));
% F1flight = zeros(size(t1flight));

%Flight Phase 2
Yapex2 = Y2(end) + V2(end)^2/(2*g);
T2apex = V2(end)/g; %time to reach apex
T2drop = abs(2*(Y01 - Yapex2)/V01); %time to drop from apex
T2flight = T2apex + T2drop;

T2 = T2 + t1flight(end);

t2flight = linspace(0,T2flight,51)' + T2(end);
x2flight = X2(end) + U2(end)*(t2flight - T2(end));
y2flight = Y2(end) - 1/2*g*(t2flight - T2(end)).^2 + V2(end)*(t2flight-T2(end));
u2flight = U2(end)*ones(size(t2flight));
v2flight = V2(end) - g*(t2flight - T2(end));
% F2flight = zeros(size(t2flight));

% combine time, states and forces into one vector each
t = [T1; t1flight(2:end); T2(2:end); t2flight(2:end)];
x = [X1; x1flight(2:end); X2(2:end) + d; x2flight(2:end) + d];
y = [Y1; y1flight(2:end); Y2(2:end); y2flight(2:end)];
u = [U1; u1flight(2:end); U2(2:end); u2flight(2:end)];
v = [V1; v1flight(2:end); V2(2:end); v2flight(2:end)];

% F1 = F1 + Fspring1;
% F2 = F2 + Fspring2;
% F = [F1; F1flight(2:end); F2(2:end); F2flight(2:end)];


n = 120;
ti = linspace(0,t(end),n+1);
xi = interp1(t,x,ti);
yi = interp1(t,y,ti);
ui = interp1(t,u,ti);
vi = interp1(t,v,ti);


%% COMPUTE ENERGY GAIN
Ef = ze2(2) + 1/2*sum(ze2(3:4).^2); % final energy
E0 = h01 + 1/2*U01.^2;
% disp(['Energy gained: ', num2str(Ef-E0)])

end