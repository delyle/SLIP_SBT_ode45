function plotSLIPsim(ti,xi,yi,L0,g,T1,T2,t1flight,t2flight,vbelt1,vbelt2)
% This function plots a freeze frame of the solved SLIP from SLIPsim.m

%% FREEZE FRAME
n = 120;
FR = 60;

% ti = linspace(0,t(end),n+1);
% xi = interp1(t,x,ti);
% yi = interp1(t,y,ti);
% ui = interp1(t,u,ti);
% vi = interp1(t,v,ti);
% Fi = interp1(t,F,ti);
li = sqrt(xi.^2 + yi.^2);

% figure('Color','white');
hold on
Lmax = L0;
plot(Lmax*[-1 1],[0 0],'k-')
yl = [0 2*Lmax];
xl = Lmax*[-5 5];

for i = 1:1:n+1
%     if i > 1
%         delete(lh)
%     end
    
    % Use magnitude of force vector to modify width of leg lines.
    alpha = min([li(i)/g,1]);
    lwidth = 4*alpha+ 0.01;
    if ti(i) <= T1(end)
        footx = 0 - vbelt1*ti(i);
        lcolor = [ti(i)/T1(end) 0 0];
        beltspeed = vbelt1;        
        lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
    elseif ti(i) >= T2(1) && ti(i) <=t2flight(1)
        footx = d - vbelt2*(ti(i) - T2(1));
        lcolor = [0 0 (ti(i)-T2(1))/(T2(end)-T2(1))];
        beltspeed = vbelt2;
        lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
    end
%     lh(1) = plot(app.UIaxes,xi(i),yi(i),'o','markersize',8,'markerfacecolor',lcolor,'markeredgecolor',lcolor);
    lh(3) = patchline([x1(1) x1(2)],[0 0],'edgecolor',[0 0 0]);
    
    xlim(xl)
    ylim(yl)
    %     set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YTick',[]);
    %     if plotvec
    %         hold on
    %         xr = xi(i) - footx;
    %         l = sqrt(xr^2 + yi(i)^2);
    %         Fx = Fi(i)*xr/l;
    %         Fy = Fi(i)*yi(i)/l;
    %
    %         lh(3) = quiver(xi(i),yi(i),ui(i)+beltspeed,vi(i),'color',lcolor,'linewidth',2,'maxheadsize',0.4);
    %         lh(4) = quiver(xi(i),yi(i),ui(i),vi(i),'color',lcolor,'linewidth',2,'linestyle','--','maxheadsize',0.4);
    %         lh(5) = quiver(xi(i),yi(i),Fx,Fy,'k','linewidth',2,'maxheadsize',0.5,'autoscalefactor',0.1);
    %     end
    drawnow
%     if saveanimation
%         % Get an image of the plot
%         frame = getframe(gcf);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if i == 1
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/FR);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/FR);
%         end
%     else
%         pause(0.03)
%     end
    
end
plot(xi,yi,'m--','LineWidth',2);
end