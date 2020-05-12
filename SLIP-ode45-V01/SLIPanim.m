function SLIPanim(ti,xi,yi,ui,vi,T1,T2,t2flight,d,vbelt1,vbelt2,saveanimation,plotvec)
%UNTITLED Summary of this function goes here
% CAREFUL: discrep is an output of SLIPsim but not a required input for
% SLIPanim.m!

%   Detailed explanation goes here

if nargin < 12
    saveanimation = true;
    plotvec = true;
end

if saveanimation
    filename = 'SLIPSBT.gif';
end


figure('Color','white');
box on
yl = [0 5];
xl = [-20 15];
n = length(ti)-1;
FR = 60;
hold on

for i = 1:n+1
    if i > 1
        delete(lh)
    end
   
    lwidth = 4;
    if ti(i) <= T1(end)
        footx = 0 - vbelt1*ti(i);
        lcolor = [1 0 0];
        beltspeed = vbelt1;
        lh(2) = plot([footx xi(i)],[0 yi(i)],'Color',lcolor,'LineWidth',lwidth);
    elseif ti(i) >= T2(1) && ti(i) <=t2flight(1)
        footx = d - vbelt2*(ti(i) - T2(1));
        lcolor = [0 0 1];
        beltspeed = vbelt2;
        lh(2) = plot([footx xi(i)],[0 yi(i)],'Color',lcolor,'LineWidth',lwidth);
    end
    lh(1) = plot(xi(i),yi(i),'o','markersize',10,'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 0]);
   
    %     alpha = min([li(i),1]);
    %     lwidth = 4*alpha+ 0.01;
    %     lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
    xlim(xl)
    ylim(yl)
    set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YTick',[]);
    if plotvec
        hold on
        xr = xi(i) - footx;
        l = sqrt(xr^2 + yi(i)^2);
       
        lh(3) = quiver(xi(i),yi(i),ui(i)+beltspeed,vi(i),'color',lcolor,'linewidth',2,'maxheadsize',1);
        lh(4) = quiver(xi(i),yi(i),ui(i),vi(i),'color',lcolor,'linewidth',2,'linestyle','--','maxheadsize',1);
    end
    drawnow
    if saveanimation
        % Get an image of the plot
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/FR);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/FR);
        end
        %     else
        %         pause(0.03)
    end
   
end

end

