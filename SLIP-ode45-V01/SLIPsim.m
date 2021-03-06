function [ti,xi,yi,ui,vi,T1,T2,t2flight,d,discrep] = SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2)
% This function solves SLIP on a split-belt treadmill
% This was the function run by SLIP_GUI to solve SLIP on a treadmill

% Previous function header
% function [ti,xi,yi,ui,vi,L0,g,T1,T2,t2flight,d,Ef,E0,xprev,yprev] = SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2)

%% Function Outputs Explained
% - ti through vi are vectors of interpolated values (ti has uniform
% intervals) of the solution
% - T1, T2, t2flight, d are for plotting purposes only
% - specifically:
% T1 is duration of Phase 1
% T2 is duration of Phase 2
% t2flight is the duration of flight phase after Phase 2
% d is the distance from the origin in Phase 1 to the origin of Phase 2.
% The origin is the point of foot contact at the beginning of each phase.

%% Possible Solution Errors

    
%% PROBLEM INFO, EQUATIONS OF MOTION
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

%% PROBLEM PARAMETERS
% Non-dimensionalization: m = g = l = 1
g = 1;
L0 = 1;
discrep = 1000;

% event parameters
Tmax = 40;
f = 1.01; %this is for failure event1 (mass falling through floor)

%% PHASE 1
X01 = L0*sin(theta01);
Y01 = L0*cos(theta01);
Z01 = [X01 Y01 U01 -V01]; % define initial conditions

% this is the event function
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt1,f)); 

% this solves the ode
[T1,Z1,te1,ze1,ie1] = ode45(@(T,Z) springmassODE(T,Z,K1,L0,vbelt1,c1),[0;Tmax],Z01,options);

% defining convenience vectors X1,Y1,U1,V1 from Z1
X1 = Z1(:,1); Y1 = Z1(:,2); U1 = Z1(:,3); V1 = Z1(:,4);

%if it hit the ground, end solution
% if ze1(end,2) < 0.1
%     n = 120; %arbitrary number
%     ti = linspace(0,T1(end),n+1);
%     xi = interp1(T1,X1,ti);
%     yi = interp1(T1,Y1,ti);
%     ui = interp1(T1,U1,ti);
%     vi = interp1(T1,V1,ti);
%     kinematicError()
%     T2 = 0;
%     t2flight = 0;
%     d = 0;
%     return
% end

%if it's too low coming out of phase 1 (based on theta02), end solution
    % not sure about this one yet since what if it's lower coming out of
    % phase 1, but has enough vertical velocity to get back up to height?
    % If you find the theoretical peak of the flight phase and it's lower,
    % then it shouldn't work. Okay... but then what? Do you take the
    % solution at the end of Phase or after the fligth phase? What does
    % this solution represent? Like, what exactly would be the use of that
    % solution? It could lead me to erroneously interpreting that solution
    % as somethign that it's not..? Not clear yet.


%% PHASE 2
X02 = L0*sin(theta02);
Y02 = L0*cos(theta02);

U02 = U1(end); % U is conserved
Yapex1 = Y1(end) + V1(end)^2/(2*g); % peak of flight phase arc
V02 = sqrt(2*g*abs(Y02 - Yapex1)); %make negative since we know it is downwards
Z02 = [X02 Y02 U02 -V02];

% this is the event function
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt2,f));

% this solves the ode
[T2,Z2,te2,ze2,ie2] = ode45(@(T,Z) springmassODE(T,Z,K2,L0,vbelt2,c2),[0;Tmax],Z02,options);
X2 = Z2(:,1); Y2 = Z2(:,2); U2 = Z2(:,3); V2 = Z2(:,4);

% failure events should have been caught from state 1 as well as not having enough height to get to theta02
% now, track same failures (hit ground, not high enough for theta01) then
% fail the solution
% Failures should be caught by error minimization function for after Phase
% 2 (or flight phase of Phase 2?)

%% JOIN PHASES INTO CONTINUOUS VECTOR

% q stands for query (made to evenly space time intervals)
Tq1 = linspace(0,T1(end),101);
Tq2 = linspace(0,T2(end),101);
% interpolate all states
Zq1 = interp1(T1,Z1,Tq1,'spline');
Zq2 = interp1(T2,Z2,Tq2,'spline');

%??? not sure if i need this anymore
%Flight phase prior to Phase 1
% Tprev = abs(V01/g);
% tprev = linspace(0,Tprev);
% xprev = X01 - U01*tprev(end) + U01.*tprev;
% yprev = h01 - 1/2*g*tprev.^2;

%Redefine time and states with queried values
%At this stage T2(1) = 0, which must be converted later
T1 = Tq1'; T2 = Tq2';
X1 = Zq1(:,1); Y1 = Zq1(:,2); U1 = Zq1(:,3); V1 = Zq1(:,4);
X2 = Zq2(:,1); Y2 = Zq2(:,2); U2 = Zq2(:,3); V2 = Zq2(:,4);

%Time to add flight phases
%Flight (after) Phase 1
T1apex = V1(end)/g; %time to reach apex (-ve cancels out with g = 1 instead of -1)
T1drop = abs(2*(Y02 - Yapex1)/abs(V02)); %time to drop from apex
T1flight = T1apex + T1drop;

d = X1(end) - X02 + U1(end)*T1flight; %"step length" to reset the origin

%convert flight-phase 1 vectors into global time
t1flight = linspace(0,T1flight,51)' + T1(end);
x1flight = X1(end) + U1(end)*(t1flight - T1(end));
y1flight = Y1(end) - 1/2*g*(t1flight - T1(end)).^2 + V1(end)*(t1flight-T1(end));
u1flight = U1(end)*ones(size(t1flight));
v1flight = V1(end) - g*(t1flight - T1(end));
% F1flight = zeros(size(t1flight));

%Flight Phase 2
Yapex2 = Y2(end) + V2(end)^2/(2*g); %apex height
T2apex = V2(end)/g; %time to reach apex
T2drop = abs(2*(Y01 - Yapex2)/V01); %time to drop from apex
T2flight = T2apex + T2drop;

%convert to flight phase 2 vectors into global time
T2 = T2 + t1flight(end);

t2flight = linspace(0,T2flight,51)' + T2(end);
x2flight = X2(end) + U2(end)*(t2flight - T2(end));
y2flight = Y2(end) - 1/2*g*(t2flight - T2(end)).^2 + V2(end)*(t2flight-T2(end));
u2flight = U2(end)*ones(size(t2flight));
v2flight = V2(end) - g*(t2flight - T2(end));
% F2flight = zeros(size(t2flight));

% combine time, states into one vector each
t = [T1; t1flight(2:end); T2(2:end); t2flight(2:end)];
x = [X1; x1flight(2:end); X2(2:end) + d; x2flight(2:end) + d];
y = [Y1; y1flight(2:end); Y2(2:end); y2flight(2:end)];
u = [U1; u1flight(2:end); U2(2:end); u2flight(2:end)];
v = [V1; v1flight(2:end); V2(2:end); v2flight(2:end)];

% CREATE INTERPOLATED VECTORS FOR FUNCTION OUTPUT
n = 120; %arbitrary number
ti = linspace(0,t(end),n+1);
xi = interp1(t,x,ti);
yi = interp1(t,y,ti);
ui = interp1(t,u,ti);
vi = interp1(t,v,ti);

kinematicError();

% Function to compute sum of absolute errors in kinematics
    function kinematicError
        discrep = sum([abs(xi(end)-xi(1)) abs(yi(end)-yi(1)) abs(ui(end)-ui(1)) abs(vi(end)-vi(1))]);
    end


end