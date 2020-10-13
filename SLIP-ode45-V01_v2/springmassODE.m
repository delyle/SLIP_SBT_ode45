function dZdT = springmassODE(T,Z,K,c)
% The set of differential equations to be submitted to ode45()
% Solves the dynamics of SLIP on a treadmill

X = Z(1);
Y = Z(2);
U = Z(3);
V = Z(4);

L = sqrt(X.^2 + Y.^2);

Ldot = (X*U+Y*V)/L;

Fspring = K*(1 - L); % L0 - L, but L0 is 1
Fdamping = -c*Ldot;

dZdT = zeros(4,1);
dZdT(1) = U;
dZdT(2) = V;
dZdT(3) = (Fspring + Fdamping).*X/L;
dZdT(4) = (Fspring + Fdamping).*Y/L - 1;