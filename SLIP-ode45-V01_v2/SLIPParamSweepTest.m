K = [5,10];
C = linspace(0,1,21);
theta0 = linspace(-pi/3,pi/3,10)+pi/2;
u_variation = linspace(-0.51,0.51,21);
vbelt1 = 0.51;
vbelt2 = 0.17;
U0 = unique(sort([vbelt1+u_variation, vbelt2+u_variation]));
V0 = linspace(-0.5,0,21);

[U0mat,V0mat,theta0mat,Kmat,Cmat] = ndgrid(U0,V0,theta0,K,C);

[Tfmat,Ufmat,Vfmat,thetafmat,flag] = SLIPParamSweep(U0mat,V0mat,theta0mat,Kmat,Cmat);

%% Get DOUBLE the data by taking advantage of the symmetry of the problem and flipping some signs

% when U0 is negative, the results are the same if theta0 is flipped in the
% y axis

% V0, Tf, Vf and flags are all the same. Uf switches sign, and thetaf is
% flipped in the y axis.


%% save results

save([date_prefix('yyyymmddHHMM'),'SweepEndpoints.mat'],'*mat','flag')

