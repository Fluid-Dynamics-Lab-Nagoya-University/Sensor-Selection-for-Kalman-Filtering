% function to select sensor based on KF
% F_sensor_KF.m : ------------- Programer : Kumi Nakai, Taku Nonomura 2020
%                               Last modified: 2021/11/17 K.Nakai(maintanence)
% Noise (system&observation): White noise is produced using amplitude given as input param
% Component number: 1

function [time, H, isensor, ranki]=F_sensor_SDP_KF(Aorg,Borg,Corg,p,sigma_s2,sigma_o2)%,Dorg

    tic
    H=[];
    [n,r]= size(Corg);
    Ctmp=Corg;
    AAI = kron(Aorg',Aorg')-eye(r^2);
    AAIinv = pinv(AAI);
    for pp=1:1:p
        if pp ~= 1
            Cprev=H*Corg;
        else
            Cprev=[];
            H=zeros(1,n);
        end
        %% Search the maximum rank of observability gramian
        maxrank=0;
        for i=1:1:n
            % Update C for next sensor selection
            C = [ Cprev ; Ctmp(i,:)];
            CTC = C'*C;
            vecCTC = reshape(CTC,[],1);
            vecWo = - AAIinv * vecCTC; % fastest calculation (200726)
            G(:,:,i) = reshape(vecWo,[r,r]);
            ranki(i) = rank(G(:,:,i)); 
            maxrank = max(rank(G(:,:,i)),maxrank);
        end
        % ===The safest approach?===
        indx = (r-maxrank+1):r;
        
        %% Sensor selection based on Kalman filter
        obj2(1:n)=10^100;
        for i=1:1:n
            if( ranki(i)==maxrank )
            %% Update C and system for next sensor selection
                C = [ Cprev ; Ctmp(i,:)];

            %% Extract observable subspace
                [ Abar, ~, Cbar, T_obsv, ~ ] = obsvf( Aorg, Borg, C ); 

                % ===the simplest approach but doesn't work as expexted===
                % for j=1:r %j=1:size(Cbar,2)
                %    Cabs(j)=norm(Cbar(:,j)); % calculate norm (for loop is needed when R2013a is used)
                % end
                % indx=find(Cabs); %element numbers whose component is nonzero are input
                    % ---> the accuracy of obsvf is not enough...
                % rs  =nnz (Cabs); %the number of nonzero elements = the rank of observability gramian
                % ===give the threshold to detect the observability===
                % indx=find(Cabs>1e-8); % ---> rounding error...
                % rs  =size(indx,1);

                Asub = Abar(indx,indx);
                Csub = Cbar( :  ,indx); 
                
            %% Noise setting
                I = eye( size (Asub) ); 
                Q = sigma_s2 * eye( r ); 
                Qtmp = T_obsv * Q * T_obsv'; 
                Qsub = Qtmp(indx,indx); 
                R = sigma_o2 * eye( pp ); 
                S = zeros( size(Csub') ); 
                
            %% Steady covariance matrix (Steady Kalman filter)
                [Pst,~,~] = idare(Asub', Csub', Qsub, R, S, I); %B*Q*B' (R2019a-)
                obj1(i) = ranki(i);
                obj2(i) = trace(Pst);
            end
        end
        [~,isensor(pp)] = min(obj2);
        H(pp,isensor(pp)) = 1;
        text = [ 'Result(KF-based): isensor=' num2str(isensor(pp)) ', rank(Wo)=' num2str(obj1(isensor(pp))) ]; disp(text);
        Ctmp(isensor(pp),:) = zeros(1,r);

    end
    time=toc;
    isensor=isensor';

end

