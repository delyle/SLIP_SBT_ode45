limitCycle = min(discrep,[],'all');

for a1 = 1:length(theta01)
    for b1 = 1:length(theta02)
        for c1 = 1:length(c)
            for d1 = 1:length(K)
                for e1 = 1:length(U01)
                    for f1 = 1:length(h01)
                        if discrep(a1,b1,c1,d1,e1,f1) == limitCycle
                            fprintf('i = %.f, j = %.f, k = %.f, l = %.f, m = %.f, n = %.f \n',a1,b1,c1,d1,e1,f1);
                            fprintf('theta01 = %f \n',theta01(a1));
                            fprintf('theta02 = %f \n',theta02(b1));
                            fprintf('c = %f \n',c(c1));
                            fprintf('K = %f \n',K(d1));
                            fprintf('U01 = %f \n',U01(e1));
                            fprintf('h01 = %f \n',h01(f1));
                            index1 = a1;
                            index2 = b1;
                            index3 = c1;
                            index4 = d1;
                            index5 = e1;
                            index6 = f1;
                        end                            
                    end
                end
            end
        end
    end
end
