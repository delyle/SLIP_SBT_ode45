function [F,Tf,Uf,Vf,thetaf] = SLIP_LookupFunctions(data,savename)
if nargin < 2
    savename = [];
end

A = load(data);

Tf = griddedInterpolant(A.U0mat,A.V0mat,A.theta0mat,A.Kmat,A.Cmat,A.Tfmat);
Uf = griddedInterpolant(A.U0mat,A.V0mat,A.theta0mat,A.Kmat,A.Cmat,A.Ufmat);
Vf = griddedInterpolant(A.U0mat,A.V0mat,A.theta0mat,A.Kmat,A.Cmat,A.Vfmat);
thetaf = griddedInterpolant(A.U0mat,A.V0mat,A.theta0mat,A.Kmat,A.Cmat,A.thetafmat);

F = @(U0,V0,theta0,K,C) [Tf(U0,V0,theta0,K,C),Uf(U0,V0,theta0,K,C),Vf(U0,V0,theta0,K,C),thetaf(U0,V0,theta0,K,C)];

if ischar(savename) && ~isempty(savename)
   save(savename,'F','Tf','Uf','Vf','thetaf')
end



