clear; close all; clc;

%Read in values
vbelt1 = 1.5;
vbelt2 = 0.5;
theta01 = linspace(-pi/3,0,21);
theta02 = linspace(-pi/3,pi/3,21);
U01 = linspace(0,5,5);
h01 = linspace(10,15,5);
% U01 = 3;
% h01 = 14;
c = linspace(0,1,21);
K = linspace(1,5,5);


discrep = zeros(length(theta01),length(theta02),length(c),length(K),length(U01),length(h01));

for i = 1:length(theta01)
    for j = 1:length(theta02)
        fprintf('i = %.f, j = %.f \n',i,j);
        for k = 1:length(c)
            for l = 1:length(K)
                for m = 1:length(U01)
                    for n = 1:length(h01)                        
                        %Solve
                        [ti,xi,yi,ui,vi,~,~,T1,T2,t2flight,d,Ef,E0,xprev,yprev] = SLIPsim(vbelt1,vbelt2,theta01(i),theta02(j),U01(m),h01(n),c(k),K(l),K(l));
                        discrep(i,j,k,l,m,n) = sum([abs(Ef-E0) abs(xi(end)-xi(1)) abs(yi(end)-yi(1)) abs(ui(end)-ui(1)) abs(vi(end)-vi(1))]);
%                         fprintf('i = %.f, j = %.f, k = %.f, l = %.f, m = %.f, n = %.f \n',i,j,k,l,m,n);
                    end
                end
            end
        end
    end
end
