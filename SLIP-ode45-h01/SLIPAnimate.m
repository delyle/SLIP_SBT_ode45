clear; close all; clc;

vbelt1 = 1.5;
vbelt2 = 0.5;
theta01 = -0.628319;
theta02 = 0.418879;
K = 5;
c = 0.15;
U01 = 2.5;
h01 = 11.25;

saveanimation = true;
plotvec = true;

if saveanimation
    filename = 'SLIPSBT.gif';
end

[ti,xi,yi,ui,vi,~,~,T1,T2,t2flight,d,Ef,E0,xprev,yprev] = SLIPsim(vbelt1,vbelt2,theta01,theta02,U01,h01,c,K,K);
n = length(ti)-1;

%% Force plots

Edotstor = zeros(length(ti),1);
Edotdisp = zeros(length(ti),1);
Edottotal = zeros(length(ti),1);


L0 = 10;
% Vtotal = sqrt(ui.^2 + vi.^2);

for i = 1:n+1
    if ti(i) <= T1(end)
        Xbelt = (xi(i) + vbelt1*ti(i)); %belt frame
        L = sqrt(Xbelt.^2 + yi(i).^2);
        Fspringx = K*(L0 - L).*Xbelt/L;
        Fspringy = K*(L0 - L).*yi(i)/L;
        Edotstor(i) = Fspringx*ui(i) + Fspringy*vi(i); %F dot product velocity
       
        Ldot = ((Xbelt*ui(i))+yi(i)*vi(i))/L;
        Fdampx = -c*Ldot.*Xbelt/L;
        Fdampy = -c*Ldot.*yi(i)/L;
        Edotdisp(i) = Fdampx*ui(i) + Fdampy*vi(i);
       
        Edottotal(i) = Edotdisp(i) + Edotstor(i);
    elseif ti(i) >= T2(1) && ti(i) <=t2flight(1)
        Xbelt = (xi(i) + vbelt1*ti(i)); %belt frame
        L = sqrt(Xbelt.^2 + yi(i).^2);
        Fspringx = K*(L0 - L).*Xbelt/L;
        Fspringy = K*(L0 - L).*yi(i)/L;
        Edotstor(i) = Fspringx*ui(i) + Fspringy*vi(i); %F dot product velocity
       
        Ldot = ((Xbelt*ui(i))+yi(i)*vi(i))/L;
        Fdampx = -c*Ldot.*Xbelt/L;
        Fdampy = -c*Ldot.*yi(i)/L;
        Edotdisp(i) = Fdampx*ui(i) + Fdampy*vi(i);
       
        Edottotal(i) = Edotdisp(i) + Edotstor(i);
    end
end

figure('Color','white');
hold on


Edisp = trapz(ti,Edotdisp);
disp(Edisp);


t1endindex = max(find(ti < T2(1))); %=71

%Truncate fast/slow leg
for i = 1:length(ti)
    if i < t1endindex
        Edottotal1(i) = Edottotal(i);
        Edotstor1(i) = Edotstor(i);
        Edotdisp1(i) = Edotdisp(i);
        ti1(i) = ti(i);
    else
        j = i - t1endindex + 1;
        Edottotal2(j) = Edottotal(i);
        Edotstor2(j) = Edotstor(i);
        Edotdisp2(j) = Edotdisp(i);
        ti2(j) = ti(i);
    end
end
plot(ti1,Edottotal1,'r','LineWidth',2);
plot(ti2,Edottotal2,'b','LineWidth',2);
plot(ti1,Edotdisp1,'r--','LineWidth',2);
plot(ti2,Edotdisp2,'b--','LineWidth',2);

%% Animate
figure('Color','white');
box on
yl = [0 20];
xl = [-35 15];
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


%% CoP Plotting
figure();
plot(xi,yi);
% f1 = -aux.vbelt1*t1 - x1;
% f2 = -aux.vbelt2*t2 + d - x2;
% f1 = f1*1000*0.88; f2 = f2*1000*0.88; %dimensionalize and convert to mm
% %Calculate SLA
% SLfast = f1(1)-f2(end);
% SLslow = f2(1)-f1(end);
% SimSLA = (SLfast - SLslow)./(SLfast + SLslow);


%Plot
% hold('on')
% n = 120;
% for i = 1:1:n+1
%     lwidth = 1;
%     if ti(i) <= T1(end)
%         footx = 0 - vbelt1*ti(i);
%         lcolor = [ti(i)/T1(end) 0 0];
%         lh(2) = plot([footx xi(i)],[0 yi(i)],'Color',lcolor,'LineWidth',lwidth);
%     elseif ti(i) >= T2(1) && ti(i) <=t2flight(1)
%         footx = d - vbelt2*(ti(i) - T2(1));
%         lcolor = [0 0 (ti(i)-T2(1))/(T2(end)-T2(1))];
%         lh(2) = plot([footx xi(i)],[0 yi(i)],'Color',lcolor,'LineWidth',lwidth);
%     end
%     lh(1) = plot(xi(i),yi(i),'o','markersize',8,'markerfacecolor',lcolor,'markeredgecolor',lcolor);
% end