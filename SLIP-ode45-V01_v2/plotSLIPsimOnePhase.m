function plotSLIPsimOnePhase(T,Z,c,K,vbelt)

% Extract variables

X = Z(:,1);
Y = Z(:,2);
U = Z(:,3);
V = Z(:,4);

% get spring and damping forces
L = sqrt(X.^2 + Y.^2);
Ldot = (X.*U+Y.*V)./L;

Fspring = K*(1 - L); % L0 - L, but L0 is 1
Fdamping = -c*Ldot;

% correct for belt speed, for the sake of kinematics in lab frame
foot = -vbelt*T;
Xlab = X + foot;
Ulab = U - vbelt;

% initialize figure
figure('color','w','position',[440   378   560   750])

%%%% Kinematics %%%%

%%% Belt Frame %%%
subplot(2,2,1)
kinematics_plot(zeros(size(X)),X,Y,U,V)

%%% Lab Frame %%%
subplot(2,2,2)
kinematics_plot(foot,Xlab,Y,Ulab,V)



% plot forces
subplot(2,2,[3 4])
plot(T,[Fspring,Fdamping])
legend('Spring Force','Damping')

end

function kinematics_plot(foot,X,Y,U,V)

plot(X,Y,'--')
hold on

% initial and final legs
plot([foot(1) X(1)],[0 Y(1)],'k-','linewidth',3) % initial
plot([foot(end) X(end)],[0 Y(end)],'k-','linewidth',3) % final

% initial and final positions
i = [1, length(X)];
plot(X(i),Y(i),'o','markersize',10) % mass

% initial and final velocities
quiver(X(i),Y(i),U(i),V(i))

axis equal
xlim([-1.1 1.1])
ylim([0 2.2])

end
