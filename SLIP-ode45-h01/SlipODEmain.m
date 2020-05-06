% This was the script used to develop the working SLIPsim.m used in SLIPGUI
% SLIPsim.m is no longer identical to this since it had to be changed to
% accomodate its use as a function, separating the solving and plotting
% into two separate functions, namely, SLIPsim.m and plotSLIPsim.m

%% INFO
% This code solves SLIP on a treadmill
% Z = [x,y,u,v];
% x,u,y,v,theta
% L(t) = sqrt(X^2 + Y^2)
% dX = U;
% dU = K*(L0 - L)*sin(theta);
% sin(theta) = x/L;
% cos(theta) = y/L;
% dY = V;
% dV = K*(1-L)*cos(theta) - 1;

close all; clc; clear;
saveanimation = true;
filename = "SLIP.gif";

%% PROBLEM DATA
% Non-dimensionalization: m = g = l = 1
g = 1;
vbelt1 = 0; % Belt speed (leftward)
vbelt2 = 0;
E0 = 20; % initial energy; constrains V0
Tmax = 40;
K = 1;
L0 = 10;
f = 1.01;

%% PHASE 1
theta01 = 0; % angle of attack (radians, measured to vertical)
X01 = L0*sin(theta01);
Y01 = L0*cos(theta01);
U01 = 0; % landing horizontal speed (rightward)
V01 = sqrt(2*(E0 - Y01) - U01^2);
%Check if Vland is valid
if ~isreal(V01)
    Emin = Y01 +  1/2*U01^2;
    txt = ['Energy must exceed ',num2str(Emin),' to allow real vertical speed'];
    error(txt)
end
Z01 = [X01 Y01 U01 -V01];
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt1,f)); % this is the event function
[T1,Z1,te1,ze1,ie1] = ode45(@(T,Z) springmassODE(T,Z,K,L0,vbelt1),[0;Tmax],Z01,options); % this solves the ode
X1 = Z1(:,1); Y1 = Z1(:,2); U1 = Z1(:,3); V1 = Z1(:,4);

%% PHASE 2
theta02 = 0;
X02 = L0*sin(theta02);
Y02 = L0*cos(theta02);
U02 = U1(end); % U is conserved
Yapex1 = Y1(end) + V1(end)^2/(2*g);
V02 = sqrt(2*g*abs(Y02 - Yapex1)); %make negative since we know it is downwards
Z02 = [X02 Y02 U02 -V02];
options = odeset('Events',@(T,Z) springmassTerminate(T,Z,L0,vbelt2,f)); % this is the event function
[T2,Z2,te2,ze2,ie2] = ode45(@(T,Z) springmassODE(T,Z,K,L0,vbelt2),[0;Tmax],Z02,options); % this solves the ode
X2 = Z2(:,1); Y2 = Z2(:,2); U2 = Z2(:,3); V2 = Z2(:,4);

%% JOIN INTO CONTINUOUS VECTOR
Tq1 = linspace(0,T1(end),101);
Tq2 = linspace(0,T2(end),101);
Zq1 = interp1(T1,Z1,Tq1,'spline');
Zq2 = interp1(T2,Z2,Tq2,'spline');

%Redefine time and states with query values
T1 = Tq1'; T2 = Tq2';
X1 = Zq1(:,1); Y1 = Zq1(:,2); U1 = Zq1(:,3); V1 = Zq1(:,4);
X2 = Zq2(:,1); Y2 = Zq2(:,2); U2 = Zq2(:,3); V2 = Zq2(:,4);

%Flight Phase 1
T1apex = V1(end)/g; %time to reach apex
T1drop = abs(2*(Y02 - Yapex1)/abs(V02)); %time to drop from apex
T1flight = T1apex + T1drop;

d = X1(end) - X02 + U1(end)*T1flight; %"step length" to reset the origin

t1flight = linspace(0,T1flight,51)' + T1(end);
x1flight = X1(end) + U1(end)*(t1flight - T1(end));
y1flight = Y1(end) - 1/2*g*(t1flight - T1(end)).^2 + V1(end)*(t1flight-T1(end));
u1flight = U1(end)*ones(size(t1flight));
v1flight = V1(end) - g*(t1flight - T1(end));
% F1flight = zeros(size(t1flight));

%Flight Phase 2
Yapex2 = Y2(end) + V2(end)^2/(2*g);
T2apex = V2(end)/g; %time to reach apex
T2drop = abs(2*(Y01 - Yapex2)/V01); %time to drop from apex
T2flight = T2apex + T2drop;

T2 = T2 + t1flight(end);

t2flight = linspace(0,T2flight,51)' + T2(end);
x2flight = X2(end) + U2(end)*(t2flight - T2(end));
y2flight = Y2(end) - 1/2*g*(t2flight - T2(end)).^2 + V2(end)*(t2flight-T2(end));
u2flight = U2(end)*ones(size(t2flight));
v2flight = V2(end) - g*(t2flight - T2(end));
% F2flight = zeros(size(t2flight));

% combine time, states and forces into one vector each
t = [T1; t1flight(2:end); T2(2:end); t2flight(2:end)];
x = [X1; x1flight(2:end); X2(2:end) + d; x2flight(2:end) + d];
y = [Y1; y1flight(2:end); Y2(2:end); y2flight(2:end)];
u = [U1; u1flight(2:end); U2(2:end); u2flight(2:end)];
v = [V1; v1flight(2:end); V2(2:end); v2flight(2:end)];

% F1 = F1 + Fspring1;
% F2 = F2 + Fspring2;
% F = [F1; F1flight(2:end); F2(2:end); F2flight(2:end)];


%% FREEZE FRAME
n = 120;
FR = 60;

ti = linspace(0,t(end),n+1);
xi = interp1(t,x,ti);
yi = interp1(t,y,ti);
ui = interp1(t,u,ti);
vi = interp1(t,v,ti);
% Fi = interp1(t,F,ti);
li = sqrt(xi.^2 + yi.^2);

figure('Color','white');
hold on
Lmax = L0;
plot(Lmax*[-1 1],[0 0],'k-')
yl = [0 2*Lmax];
xl = Lmax*[-6 6];

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
    lh(1) = plot(xi(i),yi(i),'o','markersize',8,'markerfacecolor',lcolor,'markeredgecolor',lcolor);
    lh(3) = patchline([xl(1) xl(2)],[0 0],'edgecolor',[0 0 0]);
    
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
    else
        pause(0.03)
    end
    
end
plot(xi,yi,'m--','LineWidth',2);


%% ANIMATE
% n = 120;
% FR = 60;
% 
% ti = linspace(0,t(end),n+1);
% xi = interp1(t,x,ti);
% yi = interp1(t,y,ti);
% ui = interp1(t,u,ti);
% vi = interp1(t,v,ti);
% % Fi = interp1(t,F,ti);
% li = sqrt(xi.^2 + yi.^2);
% 
% figure('Color','white');
% Lmax = L0;
% plot(Lmax*[-1 1],[0 0],'k-')
% yl = [-2*Lmax 2*Lmax];
% xl = Lmax*[-5 3];
 
% for i = 1:n+1
%     if i > 1
%         delete(lh)
%     end
%     
%     lh(1) = plot(xi(i),yi(i),'bo','markersize',10,'markerfacecolor','b');
%     % Use magnitude of force vector to modify width of leg lines.
%     alpha = min([li(i)/g,1]);
%     lwidth = 4*alpha+ 0.01;
%     if ti(i) <= T1(end)
%         footx = 0 - vbelt1*ti(i);
%         lcolor = [1 0 0];
%         beltspeed = vbelt1;
%         lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
%     elseif ti(i) >= T2(1) && ti(i) <=t2flight(1)
%         footx = d - vbelt2*(ti(i) - T2(1));
%         lcolor = [0 0 1];
%         beltspeed = vbelt2;
%         lh(2) = patchline([footx xi(i)],[0 yi(i)],'edgecolor',lcolor,'linewidth',lwidth,'edgealpha',alpha);
%     end
%     lh(3) = patchline([xl(1) xl(2)],[0 0],'edgecolor',[0 0 0]);
%     
%     xlim(xl)
%     ylim(yl)
%     %     set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YTick',[]);
%     %     if plotvec
%     %         hold on
%     %         xr = xi(i) - footx;
%     %         l = sqrt(xr^2 + yi(i)^2);
%     %         Fx = Fi(i)*xr/l;
%     %         Fy = Fi(i)*yi(i)/l;
%     %
%     %         lh(3) = quiver(xi(i),yi(i),ui(i)+beltspeed,vi(i),'color',lcolor,'linewidth',2,'maxheadsize',0.4);
%     %         lh(4) = quiver(xi(i),yi(i),ui(i),vi(i),'color',lcolor,'linewidth',2,'linestyle','--','maxheadsize',0.4);
%     %         lh(5) = quiver(xi(i),yi(i),Fx,Fy,'k','linewidth',2,'maxheadsize',0.5,'autoscalefactor',0.1);
%     %     end
%     drawnow
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
%     
% end



%% COMPUTE ENERGY GAIN
Ef = ze2(2) + 1/2*sum(ze2(3:4).^2); % final energy
disp(['Energy gained: ', num2str(Ef-E0)])

%% PLOT
% plot(T,Z)
% legend('X','Y','U','V')

%% Original ODE ANIMATE
% figure
% Tq1 = linspace(0,T1(end),101); %time query values
% Zq1 = interp1(T1,Z1,Tq1,'spline'); %state query values
% xmin = min(Zq1(:,1)); %smallest x value
% xl = xmin*f*1.1+[-2*L0 2*L0]; %setting rough x-axis limits
% lins = linspace(0.01,1,40);
% x = L0*[-sqrt(fliplr(lins)),0,sqrt(lins)]; %vector that goes from -10 to 10 to compute y-values of arc
% y = sqrt(L0^2 - x.^2); %y values of arc
% for i = 1:length(Tq1)
%     if i > 1
%         cla
%     end
%     plot([-vbelt1*Tq1(i) Zq1(i,1)],[0 Zq1(i,2)],'b-','linewidth',2) %leg
%     hold on
%     plot(Zq1(i,1),Zq1(i,2),'ko','markersize',15','markerfacecolor','k') %CoM
%     plot(x-vbelt1*Tq1(i),y,'k-') %arc
%     xlim(xl)
%     ylim([0 2*L0])
%     pause(0.02)
% end

% %Phase 2
% Tq1 = linspace(0,T1(end),101); %time query values
% Zq1 = interp1(T1,Z1,Tq1,'spline'); %state query values
% xmin = min(Zq1(:,1)); %smallest x value
% xl = xmin*f*1.1+[-2*L0 2*L0]; %setting rough x-axis limits
% lins = linspace(0.01,1,40);
% x = L0*[-sqrt(fliplr(lins)),0,sqrt(lins)]; %vector that goes from -10 to 10 to compute y-values of arc
% y = sqrt(L0^2 - x.^2); %y values of arc
% for i = 1:length(Tq1)
%     if i > 1
%         cla
%     end
%     plot([-vbelt1*Tq1(i) Zq1(i,1)],[0 Zq1(i,2)],'b-','linewidth',2) %leg
%     hold on
%     plot(Zq1(i,1),Zq1(i,2),'ko','markersize',15','markerfacecolor','k') %CoM
%     plot(x-vbelt1*Tq1(i),y,'k-') %arc
%     xlim(xl)
%     ylim([0 2*L0])
%     pause(0.02)
% end