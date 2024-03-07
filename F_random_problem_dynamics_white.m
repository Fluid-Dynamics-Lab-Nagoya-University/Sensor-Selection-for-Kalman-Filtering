function [Aorg, Borg, Corg] = F_random_problem_dynamics_white(r,n)

% /// Random dynamical system ///            
            B_sys = 0;
            cond_AAI = 1e10;
            CNT2=0;
            while (B_sys == 0) || (cond_AAI > 1e3)
                CNT2 = CNT2+1;
                sysorg = drss(r,n); % generate random discrete test model
                % (r-th order model with one input and n outputs)
                B_sys = isstable(sysorg); % Boolean value of stability
                [Aorg,Borg,Corg,Dorg] = ssdata(sysorg);
                AAI = kron(Aorg',Aorg')-eye(r^2);
                cond_AAI = cond(AAI);
                text = [ 'CNT2=', num2str(CNT2), ', B_sys=', num2str(B_sys),...
                    ', cond(AAI)=', num2str(cond_AAI) ];
                disp(text);
            end
