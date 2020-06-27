
% Animate minimum discrep for each zeta from a parameter sweep
clc; close all; clear;

datafile = 'sweepJun2220.mat'; 

A = load(datafile);

% Extract matrices and vectors from data.
disM = A.discrepMatrix;
theta01vec = A.theta01;
theta02vec = A.theta02;
U01vec = A.U01;
V01vec = A.V01;
cvec = A.c;

% Extract other parameters
vbelt1 = A.vbelt1;
vbelt2 = A.vbelt2;
K1 = A.K1;
K2 = A.K2;
zeta = A.zeta;

% initialize vectors
sz = size(cvec);
[theta01best,theta02best,U01best,V01best] = deal(NaN(sz));


for i = 1:size(disM,5)
    
   M = disM(:,:,:,:,i);
   [discrep,j] = min(M(:).^2); % I squared the discrepancy because sometimes it is negative!!
   disp(['Discrepancy is ', num2str(discrep)])
   if ~isnan(discrep)
   [I,J,K,L] = ind2sub(size(M),j); % convert linear index to subscripts
   
   theta01 = theta01vec(I);
   theta02 = theta02vec(J);
   U01 = U01vec(K);
   V01 = V01vec(L);
   c = cvec(i);
   c1 = c; c2 = c;
   
   filename = ['zeta',num2str(zeta(i)),'_mindiscrep.gif']; % name for saving the animation
   
   [ti,xi,yi,ui,vi,T1,T2,t2flight,d,~] = ...
       SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,V01,c1,c2,K1,K2);
   SLIPanim(ti,xi,yi,ui,vi,T1,T2,t2flight,d,vbelt1,vbelt2,filename,true);
   
   % save data for plotting
   theta01best(i) = theta01;
   theta02best(i) = theta02;
   U01best(i) = U01;
   V01best(i) = V01;
   end
end
%% plot best parameters against c

close all;
figure('color','w')

subplot(2,1,1)
plot(zeta,[theta01best;theta02best])
legend('\theta_1','\theta_2')
xlabel('\zeta')
ylabel('Angle of attack [rad]')
ylim([-pi/3-0.01,pi/3+0.01])
xlim(minmax(zeta))
set(gca,'ytick',pi/6*(-2:1:2),'yticklabels',{'-\pi/3','\pi/6','0','\pi/6','\pi/3'})


subplot(2,1,2)
plot(zeta,[U01best;V01best])

legend('U_1','V_1')
xlabel('\zeta')
ylabel('Contact speed (phase 1)')
xlim(minmax(zeta))
ylim(minmax([U01vec,V01vec]))