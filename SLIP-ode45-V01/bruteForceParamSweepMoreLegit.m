clear; close all; clc;

%Read in values
vbelt1 = 1.5/sqrt(9.81*0.88);
vbelt2 = 0.5/sqrt(9.81*0.88);
theta01 = linspace(-pi/3,0,21);
theta02 = linspace(-pi/3,pi/3,21);
U01 = linspace(0,5,5);
V01 = linspace(10,15,5);
% U01 = 3;
% h01 = 14;
c1 = linspace(0,1,21);
c2 = linspace(0,1,21);
K1 = linspace(1,5,5);
K2 = linspace(1,5,5);


discrep = zeros(length(theta01),length(theta02),length(U01),length(V01),length(c1),length(c2),length(K1),length(K2));

for i = 1:length(theta01)
    for j = 1:length(theta02)
        fprintf('i = %.f, j = %.f \n',i,j);
        for k = 1:length(U01)
            for l = 1:length(V01)
                for m = 1:length(c1)
                    for n = 1:length(c2)
                        for o = 1:length(K1)
                            for p = 1:length(K2)
                                %Solve SLIP
                                [ti,xi,yi,ui,vi,T1,T2,t2flight,d,discrep] = SLIPsim(vbelt1,vbelt2,theta01(i),theta02(j),U01(k),V01(l),c1(m),c2(n),K1(o),K2(p));
                                %Compute error in kinematics e.g. |x0-xf|
%                                 discrep(i,j,k,l,m,n,o,p) = sum([abs(Ef-E0) abs(xi(end)-xi(1)) abs(yi(end)-yi(1)) abs(ui(end)-ui(1)) abs(vi(end)-vi(1))]);
                                %fprintf('i = %.f, j = %.f, k = %.f, l = %.f, m = %.f, n = %.f \n',i,j,k,l,m,n);
                            end
                        end                        
                    end
                end
            end
        end
    end
end
