clc;close all;clear

load fineLegitBruteForceOct2119.mat

c_find = 5;
K_find = 5;

ic = c == 1;
iK = K == 5;

discrep2 = discrep(:,:,ic,iK,:,:); % isolate cases for particular values of c and K
discrep2 = permute(discrep2,[1,2,5,6,3,4]); % remove excess (length==1) dimensions

% discrep2 contains discrepancies for c == c_find and K == K_find as a
% function of the following control variables:
%   1. theta01
%   2. theta02
%   3. u01
%   4. h01

% Plot the discrepancy as a function of two of the controls, while holding
% the other controls constant. For example, theta01 and theta02:

[theta01_grid, theta02_grid] = meshgrid(theta01,theta02);

ih01 = h01 == 12.5; % alternatively, specify the index directly, e.g. ih01 = 3. Play around with these values.
iU01 = U01 == 2.5;

discrepGrid = discrep2(:,:,iU01,ih01);

close all;
figure;
contourf(theta01,theta02,discrepGrid,64,'linecolor','none')
caxis([0 10])
colormap('bone')
colorbar('northoutside')
ylabel('\theta_2 [rad]')
xlabel('\theta_1 [rad]')
