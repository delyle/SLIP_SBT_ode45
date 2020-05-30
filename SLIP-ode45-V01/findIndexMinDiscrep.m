limitCycle = min(abs(discrepMatrix),[],'all');
currentZeta = 0;

for a1 = 1:length(theta01)
    for b1 = 1:length(theta02)
        for c1 = 1:length(U01)
            for d1 = 1:length(V01)
                for e1 = 1:length(zeta)                    
                        if abs(discrepMatrix(a1,b1,c1,d1,e1)) == limitCycle
                            fprintf('discrep = %f \n',discrepMatrix(a1,b1,c1,d1,e1));
                            fprintf('i = %.f, j = %.f, k = %.f, l = %.f, m = %.f\n',a1,b1,c1,d1,e1);
                            fprintf('theta01 = %f \n',theta01(a1));
                            fprintf('theta02 = %f \n',theta02(b1));
                            fprintf('U01 = %f \n',U01(c1));
                            fprintf('V01 = %f \n',V01(d1));
                            fprintf('zeta = %f \n',zeta(e1)); 
                            currentZeta = zeta(e1);
                            index1 = a1;
                            index2 = b1;
                            index3 = c1;
                            index4 = d1;
                            index5 = e1;
%                             index6 = f1;
                        end                            
                   
                end
            end
        end
    end
end
