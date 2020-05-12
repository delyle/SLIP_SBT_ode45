function plotSLIPsim(ti,xi,yi,T1,T2,t2flight,vbelt1,vbelt2,d)
% ABORTED May 12, 2020. Use SLIPanim.m now

% This function plots a freeze frame of the solved SLIP from SLIPsim.m

%% NON-DIM PARAMETERS
% have to manually change these variables across plotSLIPsim.m and
% SLIPsim.m
g = 1;
L0 = 1;

%% FREEZE FRAME
n = length(ti);

li = sqrt(xi.^2 + yi.^2);

% figure('Color','white');
hold on
Lmax = L0;
% plot(Lmax*[-1 1],[0 0],'k-')
yl = [0 2*Lmax];
xl = Lmax*[-5 5];

% for i = 1:1:n+1    
%     % Use magnitude of force vector to modify width of leg lines.
% %     alpha = min([li(i)/g,1]);
% %     lwidth = 4*alpha+ 0.01;
% 
%     if ti(i) <= T1(end)
%         footx = 0 - vbelt1*ti(i);
%         lcolor = [ti(i)/T1(end) 0 0];
%         beltspeed = vbelt1;        
%         lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
%     elseif ti(i) >= T2(1) && ti(i) <=t2flight(1)
%         footx = d - vbelt2*(ti(i) - T2(1));
%         lcolor = [0 0 (ti(i)-T2(1))/(T2(end)-T2(1))];
%         beltspeed = vbelt2;
%         lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
%     end
% %     lh(1) = plot(app.UIaxes,xi(i),yi(i),'o','markersize',8,'markerfacecolor',lcolor,'markeredgecolor',lcolor);
% %     lh(3) = patchline([x1(1) x1(2)],[0 0],'edgecolor',[0 0 0]);
%     xlim(xl)
%     ylim(yl)    
% end

plot(xi,yi,'b-');
end