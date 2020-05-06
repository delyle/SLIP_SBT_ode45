% Defines events that terminate the ODE solver in springmassODE.m
% Event 1: Touches 'arc'
% Event 2: Touches the ground

function [c, isterminal,direction] = springmassTerminate(T,Z,L0,vbelt,f)

if f < 1
    error('f must be greater than 1')
end

% First condition: Leg length is exceeded
% Second condition: mass passes through ground
c = [(Z(1) + vbelt*T).^2 + Z(2).^2 - (L0*f).^2; 
     Z(2)];

isterminal = [1;1];
direction = [0;0];
end