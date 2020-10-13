function dZdT = springmassODE(T,Z,K,L0,vbelt,c)
% The set of differential equations to be submitted to ode45()
% Solves the dynamics of SLIP on a treadmill

if nargin < 6
    c = 0;
end

X = Z(1);
Y = Z(2);
U = Z(3);
V = Z(4);

Xbelt = (X + vbelt*T); %belt frame
L = sqrt(Xbelt.^2 + Y.^2);

Ldot = ((Xbelt*U)+Y*V)/L;

Fspring = K*(L0 - L);
Fdamping = -c*Ldot;

dZdT = zeros(4,1);
dZdT(1) = U;
dZdT(2) = V;
dZdT(3) = (Fspring + Fdamping).*Xbelt/L;
dZdT(4) = (Fspring + Fdamping).*Y/L - 1;