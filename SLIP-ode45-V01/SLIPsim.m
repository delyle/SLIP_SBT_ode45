function [ti,xi,yi,ui,vi,T1,T2,t2flight,d,discrep] = SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2)
% This function solves SLIP on a split-belt treadmill
% This was the function run by SLIP_GUI to solve SLIP on a treadmill

% Previous function header
% function [ti,xi,yi,ui,vi,L0,g,T1,T2,t2flight,d,Ef,E0,xprev,yprev] = SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2)

%% To-Do

% add codes to identify error/invalidation causes
% consider maximum height

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
% IF the solution is invalid, discrep is set to NaN

% outputs initialized here so that no error occurs when invalidating
% mid-function
ti = NaN;
xi = NaN;
yi = NaN;
ui = NaN;
vi = NaN;
T1 = NaN;
T2 = NaN;
t2flight = NaN;
d = NaN;
discrep = 0;
    
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
g = 1; % non-dimensionalized
L0 = 1;
discrep = 1000;

% ode45 event parameters
Tmax = 40;
f = 1.01; % this is for failure event1 (exceed leg length) in springmassTerminate.m

%% PHASE 1
X01 = -L0*sin(theta01); %theta is positive CCW
Y01 = L0*cos(theta01);
Z01 = [X01 Y01 U01 -V01]; % define initial conditions

% this is the event function
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt1,f)); 

% this solves the ode system
[T1,Z1,te1,ze1,ie1] = ode45(@(T,Z) springmassODE(T,Z,K1,L0,vbelt1,c1),[0;Tmax],Z01,options);

% defining convenience vectors X1,Y1,U1,V1 from Z1
X1 = Z1(:,1); Y1 = Z1(:,2); U1 = Z1(:,3); V1 = Z1(:,4);

% Error catching for Phase 1
% hits ground during Phase 1
if ze1(2) <= 0 % y = 0
    invalidate();
    warning("hit ground in Phase 1")
    return
end
% peak height not sufficient to enter Phase 2
Ypeak1 = Y1(end) + (V1(end))^2/(2*g); %fill in this equation
if Ypeak1 < L0*cos(theta02) % Y02
    invalidate();
    warning("peak height too low to enter Phase 2")
    return
end

%% PHASE 2
X02 = -L0*sin(theta02);
Y02 = L0*cos(theta02);

U02 = U1(end); % U is conserved
Yapex1 = Y1(end) + V1(end)^2/(2*g); % peak of flight phase arc after Phase 1

if Yapex1 < Y02
    invalidate();
    warning("Yapex1 too low to enter Phase 2");
    return
end
if V1(end) < 0
    invalidate();
    warning("No peak in flight out of Phase 1");
    return
end

V02 = -sqrt(2*g*(Yapex1 - Y02)); % make negative since we know it is downwards
Z02 = [X02 Y02 U02 V02];

% this is the event function
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt2,f));

% this solves the ode
[T2,Z2,te2,ze2,ie2] = ode45(@(T,Z) springmassODE(T,Z,K2,L0,vbelt2,c2),[0;Tmax],Z02,options);
X2 = Z2(:,1); Y2 = Z2(:,2); U2 = Z2(:,3); V2 = Z2(:,4);

% Error catching for Phase 2
% hits ground during Phase 2
if ze2(2) <= 0 % y = 0
    invalidate();
    warning("hit ground in Phase 2")
    return
end

%% JOIN PHASES INTO CONTINUOUS VECTOR

% q stands for query (made to evenly space time intervals)
Tq1 = linspace(0,T1(end),101);
Tq2 = linspace(0,T2(end),101);
% interpolate all states for Phase 1 & 2 (no flight phases so far)
Zq1 = interp1(T1,Z1,Tq1,'spline');
Zq2 = interp1(T2,Z2,Tq2,'spline');

%Redefine time and states with queried values
%At this stage T2(1) = 0, which must be converted later
T1 = Tq1'; T2 = Tq2';
X1 = Zq1(:,1); Y1 = Zq1(:,2); U1 = Zq1(:,3); V1 = Zq1(:,4);
X2 = Zq2(:,1); Y2 = Zq2(:,2); U2 = Zq2(:,3); V2 = Zq2(:,4);

%Now add flight phases and stitch them together with Phase 1 & 2
%Flight (after) Phase 1

% only proceed if V1(end) is UPWARDS
if V1(end) < 0
    invalidate();
    warning('V1(end) downwards!');
    return
end

T1apex = V1(end)/g; %time to reach apex (-ve cancels out with g = 1 instead of -1)
T1drop = roots([1/2*g 0 Y02-Yapex1]); %time to drop from apex
T1drop = T1drop(T1drop > 0); % select positive root
if length(T1drop) > 1 % checks that there is only 1 value for T1drop
    warning('two positive roots for T1drop');
end

T1flight = T1apex + T1drop; % compute total flight time from Phase 1

d = X1(end) - X02 + U1(end)*T1flight;  % "step length" to reset the origin to the global frame

%convert flight-phase 1 vectors into global time
t1flight = linspace(0,T1flight,51)' + T1(end);
x1flight = X1(end) + U1(end)*(t1flight - T1(end));
y1flight = Y1(end) - 1/2*g*(t1flight - T1(end)).^2 + V1(end)*(t1flight-T1(end));
u1flight = U1(end)*ones(size(t1flight));
v1flight = V1(end) - g*(t1flight - T1(end));

% Flight (after) Phase 2

VdropEntry2 = 0; % Drop velocity from max height of flight Phase 2 (downwards is +ve)

if V2(end) > 0 % if V2(end) upwards
    % if there will be an APEX in Phase 2 flight
    Yapex2 = Y2(end) + V2(end)^2/(2*g); %apex height
    YflightMax2 = Yapex2; % max height during flight phase is Yapex2 if there is an apex
    T2apex = V2(end)/g; %time to reach apex

    % if apex lower than Y01, simulate flight until CoM hits ground
    YnextPhase = Y01; % target drop height from Yapex2
    if Yapex2 < Y01
        YnextPhase = 0; % just go hit the ground if can't reach Y01
    end
else    
    % there will be NO APEX in flight Phase 2
    YflightMax2 = Y2(end); % max height will be at Phase 2 exit
    T2apex = 0; % no time required to reach apex
    VdropEntry2 = -V2(end); % -ve added since downward is positive for the T2drop calcluation    
    YnextPhase = 0; %i.e. go straight to the ground
end

T2drop = roots([1/2*g VdropEntry2 YnextPhase-YflightMax2]); %time to drop from apex to YnextPhase (downwards is +ve)
% Check for complex roots and invalidate if found
if ~isreal(T2drop)
    invalidate();
    warning('T2drop is complex. This is most likely a problem with the model');
    return
end

T2drop = T2drop(T2drop > 0); % select positive root
if length(T2drop) > 1 % checks that there is only 1 value for T1drop
    warning('two positive roots for T2drop');
end

T2flight = T2apex + T2drop;

%convert to flight phase 2 vectors into global time
T2 = T2 + t1flight(end);

t2flight = linspace(0,T2flight,51)' + T2(end);
x2flight = X2(end) + U2(end)*(t2flight - T2(end));
y2flight = Y2(end) - 1/2*g*(t2flight - T2(end)).^2 + V2(end)*(t2flight-T2(end));
u2flight = U2(end)*ones(size(t2flight));
v2flight = V2(end) - g*(t2flight - T2(end));

% combine time, states into one vector each
t = [T1; t1flight(2:end); T2(2:end); t2flight(2:end)];
x = [X1; x1flight(2:end); X2(2:end) + d; x2flight(2:end) + d];
y = [Y1; y1flight(2:end); Y2(2:end); y2flight(2:end)];
u = [U1; u1flight(2:end); U2(2:end); u2flight(2:end)];
v = [V1; v1flight(2:end); V2(2:end); v2flight(2:end)];

% CREATE INTERPOLATED VECTORS FOR FUNCTION OUTPUT
n = 60; %arbitrary number
ti = linspace(0,t(end),n+1);
xi = interp1(t,x,ti);
yi = interp1(t,y,ti);
ui = interp1(t,u,ti);
vi = interp1(t,v,ti);


%% Track discrepancy @ fall point through ground (DON'T FORGET: ALSO REQUIRE STORAGE OF ZF AND ZI)

%Phase 1
% note: Tgnd1 = time to free fall to the ground from Phase 1 start
Tgnd1 = roots([g/2 V01 -Y01]); %time to fall to ground
Tgnd1 = Tgnd1(Tgnd1 >= 0); %select positive root only
Xgnd1 = U01*Tgnd1 + X01; %horz distance travelled during that time

%Phase 2
% note: Tgnd2 = time to free fall to the ground from Phase 2 end
Tgnd2 = roots([-g/2 V2(end) Y2(end)]); %time to fall to ground from beginning of flight phase 2
Tgnd2 = Tgnd2(Tgnd2 >= 0); %select positive root only
Xgnd2 = U2(end)*Tgnd2 + (X2(end) + d); %horz distance travelled during that time

discrep = Xgnd2 - Xgnd1;

%% Function to invalidate a solution
    function invalidate
        discrep = NaN;
    end

end