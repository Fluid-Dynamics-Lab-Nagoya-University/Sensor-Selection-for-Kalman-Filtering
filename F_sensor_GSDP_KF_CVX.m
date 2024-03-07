% function to select sensor based on KF
% F_sensor_KF.m : ------------- Programer : Kumi Nakai, Taku Nonomura 2020
%                               Last modified: 2021/11/17 K.Nakai(maintanence)
% Noise (system&observation): White noise is produced using amplitude given as input param
% Component number: 1

function  [zhat,z,X,Y,W,sensors]=F_sensor_GSDP_KF_CVX(A,C,Q,sigma_s2,p,lambda)
    sigma_s=sqrt(sigma_s2);
    H=sigma_s*C;
    [n,r]=size(C);
    R =eye(n)*sigma_s2;
    Rh=eye(n)*sigma_s;
    X=eye(r);
    Y=zeros(r,n);
    W=zeros(r,r);
    I=eye(r);

% cvx_begin 
%      variable z(n)
%      minimize( trace_inv(C'*diag(z)*C) )
%      ones(1,n)*z== p
%      z >= 0
%      z <= 1
% cvx_end    
%lambda=0;
%[ X-I A'-C'*Y'; A-Y*C X]
solver_now = cvx_solver;
cvx_solver('Mosek')
 cvx_begin sdp
    variable X(r,r) symmetric   
    variable Y(r,n)
    variable W(n,n) symmetric
    variable gma
    %Ynorm=vecnorm(Y,2)
    for i=1:n
        Ynorm(i,1)=norm(Y(:,i),2);
    end
    minimize( gma + lambda*ones(1,n)*Ynorm)
subject to
    - trace( Q*X ) - trace(W) + gma >= 0  %replace >= for cvx?
    [ X-I A'*X-C'*Y'; X*A-Y*C X] >= 0
    X>=0
    [ W Rh*Y'; Y*Rh X] >= 0    
cvx_end       
Ynorm=[];
Ynorm=vecnorm(Y,2,1);
z=Ynorm';
[Ynormsort,sensors]=sort(Ynorm,'descend');
thres=Ynormsort(n-p); zhat=(z>thres);
end

