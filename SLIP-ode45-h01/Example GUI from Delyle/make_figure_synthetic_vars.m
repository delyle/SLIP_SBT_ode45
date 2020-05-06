function make_figure_synthetic_vars(Forces,Times,save_name)
% User-defined values:

% % Vertical Forces
% FHRz = 7/5; % The ratio of front to hind peak vertical forces
% sh21Rzf = -0.1; % Forelimbs: The relative amplitude of the second to first harmonic
% sh31Rzf = 0.25; % Forelimbs: The relative amplitude of the third to first harmonic
% sh21Rzh = 0.1; % Hindlimbs: The relative amplitude of the second to first harmonic
% sh31Rzh = 0.5; % Hindlimbs: The relative amplitude of the third to first harmonic
% 
% % Anterior forces
% FHRy = 8/5; % The ratio of front to hind peak anterior forces
% sh42Ryf = -0.4; % Forelimbs: The relative amplitude of the fourth to second harmonic
% sh42Ryh = -0.5; % Hindlimbs: The relative amplitude of the fourth to second harmonic
% 
% Tsf = 0.65; % The duty factor of the forelimbs
% Tsh = 0.68; % The duty factor of the hindlimbs
% Ttdh = 0.28; % The touchdown time of the hindlimbs, relative to previous forelimb touchdown
%              % normalized by stride period

FHRz = Forces.Vertical.FHRz;
sh21Rzf = Forces.Vertical.sh21Rzf;
sh31Rzf = Forces.Vertical.sh31Rzf;
sh21Rzh = Forces.Vertical.sh21Rzh;
sh31Rzh = Forces.Vertical.sh31Rzh;

FHRy = Forces.Horizontal.FHRy;
sh42Ryf = Forces.Horizontal.sh42Ryf;
sh42Ryh = Forces.Horizontal.sh42Ryh;

Tsf = Times.Tsf;
Tsh = Times.Tsh;
Ttdh = Times.Ttdh;

% --- Establish forces during stance --- %:

t = (0:0.001:1)';
Fz_f = sin(pi*t)+sh21Rzf*sin(2*pi*t)+sh31Rzf*sin(3*pi*t);
Fz_h = (sin(pi*t)+sh21Rzh*sin(2*pi*t)+sh31Rzh*sin(3*pi*t))/FHRz;

Fy_f = 0.15*(-sin(2*pi*t)+sh42Ryf*sin(4*pi*t));
Fy_h = 0.15*(-sin(2*pi*t)+sh42Ryh*sin(4*pi*t))/FHRy;


% --- Translate forces into normalized stride time --- %

% First need to convert forces to the proper format:
Fz_f = [t, Fz_f]; Fz_h = [t, Fz_h];
Fy_f = [t, Fy_f]; Fy_h = [t, Fy_h];

% convert to stride time. Note that RF touchdown is t = 0;

t0R = 0; tfR = Tsf;
Fz_Rf = stance2gait(Fz_f,t0R,tfR,t); Fy_Rf = stance2gait(Fy_f,t0R,tfR,t);
t0L = 0.5; tfL = 0.5+Tsf;
Fz_Lf = stance2gait(Fz_f,t0L,tfL,t); Fy_Lf = stance2gait(Fy_f,t0L,tfL,t);
t0L = Ttdh; tfL = Ttdh+Tsh;
Fz_Lh = stance2gait(Fz_h,t0L,tfL,t); Fy_Lh = stance2gait(Fy_h,t0L,tfL,t);
t0R = 0.5 + Ttdh; tfR = 0.5+Ttdh+Tsh;
Fz_Rh = stance2gait(Fz_h,t0R,tfR,t); Fy_Rh = stance2gait(Fy_h,t0R,tfR,t);


% --- Determine total forces --- %

Fz_tot = Fz_Rf+Fz_Rh+Fz_Lf+Fz_Lh;  Fy_tot = Fy_Rf+Fy_Rh+Fy_Lf+Fy_Lh;

% scale the vertical forces so that the animal is supported
scale = mean(Fz_tot);
Fz_tot = Fz_tot/scale;
Fz_Lf = Fz_Lf/scale; Fz_Rf = Fz_Rf/scale;
Fz_Lh = Fz_Lh/scale; Fz_Rh = Fz_Rh/scale;


% --- Determine the velocities --- %

% first, determine initial velocities

Fydt = NaN(size(Fy_tot));
Fydt(1) = 0;
for i = 2:length(Fy_tot)
    Fydt(i) = trapz(t(1:i),Fy_tot(1:i));
end
vy0 = -trapz(t,Fydt);

Fzdt = NaN(size(Fz_tot));
Fzdt(1) = 0;
for i = 2:length(Fz_tot)
    Fzdt(i) = trapz(t(1:i),Fz_tot(1:i));
end
vz0 = 0 + 1/2 - trapz(t,Fzdt);

% Now calculate the CoM velocity

vy = NaN(size(t));
vy(1) = vy0;
for i = 2:length(vy)
   vy(i) = vy0 + trapz(t(1:i),Fy_tot(1:i)); 
end

vz = NaN(size(t));
vz(1) = vz0;
for i = 2:length(vz)
   vz(i) = vz0 + trapz(t(1:i),Fz_tot(1:i)) - t(i);
end

% ----- Make a figure ----- %
close(figure(1))
h = figure(1);
set(h,'position', [0,0,1500,840], 'color', 'w')


% Plot forces
subplot(4,2,1)
    hold on
    plot(t,Fz_Lf)%plot(t,Fz_Lf)
    plot(t,Fz_Rf)
    plot(t,Fz_Lh, '--')
    plot(t,Fz_Rh, '--')
    ylabel('Vertical F^*')
    box on
    legend('LF','RF','LH','RH','location', 'southwest')
    
 % total
 
subplot(4,2,3)
    hold on
    plot(t,Fz_tot)
    plot([t(1) t(end)],[1 1],'k--')
    ylabel('Total Vertical F^*')
    box on

% Plot anterior forces 
 % from each leg
subplot(4,2,5)
    hold on
    plot(t,Fy_Lf)
    plot(t,Fy_Rf)
    plot(t,Fy_Lh, '--')
    plot(t,Fy_Rh, '--')
    ylabel('Anterior F^*')
    box on
    
 % total   
subplot(4,2,7)
    hold on
    plot(t,Fy_tot)
    plot([t(1) t(end)],[0 0],'k--')
    ylabel('Total Anterior F^*')
    xlabel('t^*')
    box on
    
% Plot hodograph

subplot(4,2,[2,4])
    hold on
    plot(vy,vz);
    h(:,1) = plot(vy(1),vz(1), 'go', 'markerfacecolor', 'g'); % start marker
    u = vy(20)-vy(1);
    v = vz(20)-vz(1);
    h(:,2) = quiver(vy(1),vz(1),u,v,'color','k','maxheadsize',1);
    box on
    xlabel('v_y^* - v^*_{yAv}')
    ylabel('v^*_z')
    legend(h,{'start','direction'},'location','southeast')

% Plot COM velocities

subplot(4,2,6)
    hold on
    plot(t,vz/max(vz))
    plot(t,vy/max(vy))
    plot([t(1) t(end)],[0 0],'k--')
    ylabel('(v^* - v_{Av}^*) / v^*_{max}')
    xlabel('t^*')
    legend('V_z','V_y','location','southeast')
    box on
    
txt = [ '\sl        Vertical Forces \rm', char(10),char(10)...
        'Front-to-hind peak force ratio:  ', num2str(FHRz), char(10),...
        '\sl{Forelimbs:} \rm', char(10),...
        '   Relative amplitude of...',char(10),...
        '       2^{nd} to 1^{st} harmonic:  ' num2str(sh21Rzf), char(10),...
        '       3^{rd} to 1^{st} harmonic:  ' num2str(sh31Rzf), char(10),...
        '\sl{Hindlimbs:} \rm', char(10),...
        '   Relative amplitude of...',char(10),...
        '       2^{nd} to 1^{st} harmonic:  ' num2str(sh21Rzh), char(10),...
        '       3^{rd} to 1^{st} harmonic:  ' num2str(sh31Rzh)];
text(-0.1,-0.5,txt, 'units', 'normalized', 'VerticalAlignment', 'top', 'fontsize', 12);

txt = [ '\sl        Anterior Forces \rm', char(10),char(10)...
        'Front-to-hind peak force ratio:  ', num2str(FHRy), char(10),...
        '\sl{Forelimbs:} \rm', char(10),...
        '   Relative amplitude of...',char(10),...
        '       4^{th} to 2^{nd} harmonic:  ' num2str(sh42Ryf),char(10),...
        '\sl{Hindlimbs:} \rm', char(10),...
        '   Relative amplitude of...',char(10),...
        '       4^{th} to 2^{nd} harmonic:  ' num2str(sh42Ryh)];
text(0.35,-0.5,txt, 'units', 'normalized', 'VerticalAlignment', 'top','fontsize', 12);

txt = [ '\sl        Stance times \rm', char(10),char(10)...
        'Duty factor of forelimbs:  ', num2str(Tsf), char(10),...
        'Duty factor of hindlimbs:  ', num2str(Tsh), char(10),...
        'Touchdown time of hindlimb,', char(10),...
        ' following forelimb touchdown:  ', num2str(Ttdh)];
text(0.8,-0.5,txt, 'units', 'normalized', 'VerticalAlignment', 'top','fontsize', 12);

txt = [save_name, '.pdf'];
this_dir = cd;
cd(homedir)
pause(0.1)
export_fig(['Quad_walking_analysis/figures/',txt]);
cd(this_dir)
end