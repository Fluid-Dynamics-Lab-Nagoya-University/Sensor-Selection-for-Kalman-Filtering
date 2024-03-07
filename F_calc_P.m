% F_calc_detP.m : ------------- Programer : Kumi Nakai 2020
%                               Last modified: 2021/11/12 K.Nakai

function [Pst]=F_calc_P(p, sensors, Aorg, Borg, Corg, sigma_s2, sigma_o2)%Dorg

    [~,r]=size(Corg);
    C = Corg(sensors,:);
%     D = H*Dorg;
    CTC = C'*C;
    vecCTC = reshape(CTC,[],1);
    % vecWo = - ( kron(Aorg',Aorg')-eye(r^2) ) \ vecCTC;
    vecWo = - pinv( kron(Aorg',Aorg')-eye(r^2) ) * vecCTC;%(2021/4/22)
    Gram = reshape(vecWo,[r,r]);
    maxrank = rank(Gram);
    indx=(r-maxrank+1):r;
    [Abar, ~, Cbar, T_obsv, ~] = obsvf(Aorg, Borg, C); 
    Asub=Abar(indx,indx);
    Csub=Cbar( :  ,indx);
    %%Noise setting
    I = eye( size ( Asub ) ) ; 
    Q = sigma_s2 * eye( r );
    Qtmp = T_obsv * Q * T_obsv';
    Qsub = Qtmp(indx,indx);
    R = sigma_o2 * eye( p ) ;
    S = zeros( size( Csub' ) );
    % sum(k_obsv);
    %% Steady covariance matrix (Steady Kalman filter)
    % [Pst,~,~] =  dare(Asub', Csub', Qsub, R, S, I); %B*Q*B'
    [Pst,~,~] = idare(Asub', Csub', Qsub, R, S, I); %B*Q*B'
   % trP = trace(Pst) ;
    
end

% function [detP]=F_calc_detP(p, H,Aorg,Borg,Corg,Dorg,sigma_s2,sigma_o2)

%     [~,r]=size(Corg);
%     C = H*Corg;
%     D = H*Dorg;
%     sys = ss(Aorg,Borg,C,D,[]);
%     Gram = gram(sys,'o');
%     maxrank = rank(Gram);
%     indx=(r-maxrank+1):r;
%     [ Abar , Bbar , Cbar , T_obsv , k_obsv ] = obsvf( Aorg , Borg , C ) ; 
%     Asub=Abar(indx,indx);
%     Csub=Cbar( :  ,indx);
%     %%Noise setting
%     I = eye( size ( Asub ) ) ; 
%     Q = sigma_s2 * eye( r );
%     Qtmp = T_obsv * Q * T_obsv';
%     Qsub = Qtmp(indx,indx);
%     R = sigma_o2 * eye( p ) ;
%     S = zeros( size( Csub' ) );
%     % sum(k_obsv);
%     %% Steady covariance matrix (Steady Kalman filter)
%     % [Pst,~,~] =  dare(Asub', Csub', Qsub, R, S, I); %B*Q*B'
%     [Pst,~,~] = idare(Asub', Csub', Qsub, R, S, I); %B*Q*B'
%     detP = det(Pst) ;
    
% end