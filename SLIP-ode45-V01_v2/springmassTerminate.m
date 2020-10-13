function [c, isterminal,direction] = springmassTerminate(T,Z,f)
% Defines events that terminate the ODE solver in springmassODE.m
% Event 1: Touches 'arc'
% Event 2: Touches the ground
% Z = [X,Y,U,V];
% f is max leg length to detect mass outside arc (should be greater than
% L0=1)

if f < 1
    error('f must be greater than 1')
end

% First condition: Leg length is exceeded
% Second condition: mass passes through ground
c = [Z(1).^2 + Z(2).^2 - f.^2; 
     Z(2)];

isterminal = [1;1];
direction = [0;0];
end