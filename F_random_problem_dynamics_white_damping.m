function [Aorg, Borg, Corg] = F_random_problem_dynamics_white_damping(r,n,f,g)

% /// Linear system of three dynamical modes ///
            % specify characteristic frequencies and growth/decay rates
            % associated with continuous-time dynamics
            assert(length(f)==length(g))
            % construct low-rank continuous-time operator (rank=k)
            k = 2*length(f);  % (to perform DMD/TDMD with correct rank, set r=k)
            A1 = [];
            for ii = 1:length(f)
                A2 = [[g(ii) 2*pi*f(ii); -2*pi*f(ii) g(ii)]];
                A1 = [A1 A2];
            end
            Alowrank = [];
            for ii = 1:length(f)
                Alowrank = blkdiag(Alowrank,A1(:,(ii-1)*2+1:2*ii));
            end
            % size(Alowrank)
            dt = 0.01; % time step size
            ii = 1;
            Aorg = expm(dt*ii*Alowrank);
            AAI = kron(Aorg',Aorg')-eye(r^2);
            sysorg = drss(r,n);
            [~,Borg,Corg,Dorg] = ssdata(sysorg);
