clear
datafile = '202010131616SweepEndpoints.mat';

F = SLIP_LookupFunctions(datafile);

% for a given input to phase 1, find a phase 2 that minimizes discrepancy

% inputs
vbelt1 = 0.5;
vbelt2 = 0.17;
K = 10;
c = 1;
U01lab = 0.1;
V01 = -0.1;
theta01 = -0.5+pi/2;

Zf = F(U01lab+vbelt1,V01,theta01,K,c);

Tf1 = Zf(1);
Uf1 = Zf(2);
Vf1 = Zf(3);
thetaf1 = Zf(4);

Yf1 = sin(thetaf1);

% parameters to search
g = 1;

Yapex1 = Yf1 + (Vf1)^2/(2*g); % parabolic apex. Not reached if Vf < 0!

if Vf1 < 0
     Ypeak1flight = Yf1;
else
     Ypeak1flight = Yapex1;
end

% find range of possible initial angles

theta02min1 = pi/4;
theta02max1 = asin(Ypeak1flight);
if theta02max1 < theta02min1
    warning('No feasible solution') 
end
theta02min2 = pi - theta02max1;
theta02max2 = 3*pi/4;

theta02 = [theta02min1:0.0175:theta02max1, theta02min2:0.0175:theta02max2];
