function [Ur,Sr, Vr, Psi, S, V, Xorg, dSdiag, Rdiag, R] = F_random_problem_corr_noise(n,m,r1,r2)
    Urand=zeros(n,m);
    Vrand=zeros(m,m);
    A1=randn(n,n);
    A2=randn(m,m);
    [Q1, ~]=qr(A1);
    [Q2, ~]=qr(A2);
    Uorg=Q1(1:n,1:m);
    Vorg=Q2';
    for j=1:m
        Sorg(j,j)=1/sqrt(double(j));
    end
    Psi=Uorg(:,1:r1);
    S=Sorg(1:r1,1:r1);
    V=Vorg(:,1:r1);
    Xorg=Uorg*Sorg*Vorg';

    [~,rr]=size(Sorg);
    Ur=Uorg(:,r1+1:r2);
    Sr=Sorg(r1+1:r2,r1+1:r2);
    Vr=Vorg(:,r1+1:r2);
    Xorglr=Uorg(:,1:r1)*Sorg(1:r1,1:r1)*Vorg(:,1:r1)';
    Rlrdiag=diag(Ur*Sr*Sr*Ur');
    dSdiag=zeros(n,1);

    for i=1:n
        tmp=0;
        for j=1:rr
            tmp=tmp+(Xorg(i,j)-Xorglr(i,j))^2;
        end
        dSdiag(i)=tmp-Rlrdiag(i);
    end
    Rdiag=Rlrdiag+dSdiag;
    R=Ur*Sr*Sr*Ur'+diag(dSdiag);

