function [Tf,Uf,Vf,thetaf,flag] = SLIPParamSweep(U0,V0,theta0,K,C)


% Initialize final state matrices

[Tf,Uf,Vf,thetaf] = deal(NaN(size(U0)));

flag = false(size(U0));
wbh = waitbar(0,'Simulating SLIP');
n = numel(Tf);
tic
for i = 1:n

% perform simulation
[T,Z,~,ze1,~] = SLIPsim(U0(i),V0(i),theta0(i),K(i),C(i));

% hits ground?
if ze1(2) <= 0 % y = 0
    flag(i) = true;
end

Tf(i) = T(end);
Uf(i) = Z(end,3);
Vf(i) = Z(end,4);

Lf = sqrt(Z(end,1).^2 + Z(end,2).^2);
thetaf(i) = acos(Z(end,1)/Lf);
if mod(i,50) == 0
    % every 50 iterations, update the waitbar
    t = toc;
    tm = t/60;
    waitbar(i/n,wbh,sprintf('Simulating SLIP, %.0f minutes remaining',(n-i)/i*tm))
end

end
close(wbh)