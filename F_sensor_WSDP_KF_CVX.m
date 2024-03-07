% function to select sensor based on KF
% F_sensor_KF.m : ------------- Programer : Kumi Nakai, Taku Nonomura 2020
%                               Last modified: 2021/11/17 K.Nakai(maintanence)
% Noise (system&observation): White noise is produced using amplitude given as input param
% Component number: 1

function  [zhat,z,X,Y,Z,sensors]=F_sensor_WSDP_KF_CVX(A,C,Q,sigma_s2,p)
    sigma_s=sqrt(sigma_s2);
    H=sigma_s*C;
    [n,r]=size(C);
    Qinv=inv(Q);
    X=eye(r);
    Y=X;
    Z=zeros(r,n);
    z=p/n*ones(n,1);
    G=diag(z);
    I=eye(r);
% cvx_begin 
%      variable z(n)
%      minimize( trace_inv(C'*diag(z)*C) )
%      ones(1,n)*z== p
%      z >= 0
%      z <= 1
% cvx_end    
    size(A);
    size(H);
    Zrr=zeros(r,r);
    Zrn=zeros(r,n);
    Znr=zeros(n,r);
    solver_now = cvx_solver;
    cvx_solver('Mosek')
cvx_begin sdp
    variable X(r,r) symmetric
    variable Y(r,r) symmetric
    variable Z(r,n)
    variable z(n)
    G=diag(z);
    minimize( trace(X) )
subject to
    [X I ; I Y] >= 0d0
    [Y             Y*A-Z*H*A  Y-Z*H      Z         ;...
     A'*Y-A'*H'*Z' Y          Zrr        Zrn       ;...
     Y-H'*Z'       Zrr        Qinv       Zrn       ;...
     Z'            Znr        Znr        G         ]>=0
    ones(1,n)*z <= p
    z >= 0
    z <= 1
cvx_end
[zsort,sensors]=sort(z,'descend');
thres=zsort(n-p); zhat=(z>thres);
end

