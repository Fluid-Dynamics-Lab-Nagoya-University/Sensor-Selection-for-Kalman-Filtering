function [zhat L ztilde Utilde,sensors, data] = F_sens_sel_KF_approxnt(A, U, Q, sigma2, k, MAXITER)
% Solves the problem
%	maximize log det (sum_{i=1}^m z_i a_i a_i^T) + kappa sum_{i=1}^m(log(z_i)+ log(1-z_i))
%	subject to sum(z) = k
%			   0 <= z_i <= 1, i=1,..., 
% variable z in R^m
% problem parameters kappa (>0), a_1, ..., a_m in R^n
%
% see paper Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% Nov 2007 Sidd
% Newton's method parameters
%MAXITER  = 30;
NT_TOL = 1e-7;
GAP = 1.005;
GAP = 1.001;
%GAP = 1.0001;
% Backtracking line search parameters
alpha = 0.01;
beta = 0.5;

[m n] = size(U);

R = sigma2*eye(m);
S = zeros( size(U') );
I = eye( size(A) );

z = ones(m,1)*(k/m); % initialize
g = zeros(m,1);
ones_m = ones(m,1);
kappa = log(GAP)*n/m; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

fprintf('\nIter.  Step_size  Newton_decr.  Objective  trinv trinv_dis time \n');
%[P,~,~]=idare(A',U'*diag(sqrt(z)),Q,R,S,I);
[~,Ssvd,Vsvd]=svd(diag(sqrt(z))*U,'econ');
[P,~,~]=idare(A',Vsvd*Ssvd,Q,sigma2*eye(n),zeros(n,n),I);
fz = trace(P) - kappa*sum(log(z) + log(1-z));

fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, trace(P));
tic
idata=0;
for iter=1:MAXITER
    %tic
    %[P,~,~]=idare(A',U'*diag(sqrt(z)),Q,R,S,I);
    [~,Ssvd,Vsvd]=svd(diag(sqrt(z))*U,'econ');
    [P,~,~]=idare(A',Vsvd*Ssvd,Q,sigma2*eye(n),zeros(n,n),I);
    Pinv        = inv(P);
    Winv        = inv(Pinv+sigma2*U'*diag(z)*U);
    AWinvPinv   = A*Winv*Pinv;
    IAWinvPinv  = inv(eye(n*n)-kron(AWinvPinv,AWinvPinv));
    tIAWinvPinv = transpose(vec(eye(n)))*IAWinvPinv;
    for i=1:m
        P1(:,:,i)=-sigma2*reshape(IAWinvPinv*vec(A*Winv*U(i,:)'*U(i,:)*Winv*A'),[n,n]);
        g(i)     =trace(P1(:,:,i)) - kappa*(1./z(i) - 1./(1-z(i)));
    end
    %toc
    %tic
    for j=1:m
        for i=j:m
            L1 =  A*Winv*(-Pinv*P1(:,:,j)*Pinv+sigma2*U(j,:)'*U(j,:))*Winv*(-Pinv*P1(:,:,i)*Pinv+sigma2*U(i,:)'*U(i,:))*Winv*A';
            L2 =  A*Winv*(-Pinv*P1(:,:,i)*Pinv+sigma2*U(i,:)'*U(i,:))*Winv*(-Pinv*P1(:,:,j)*Pinv+sigma2*U(j,:)'*U(j,:))*Winv*A';
            L4 = -A*Winv*Pinv*P1(:,:,j)*Pinv*P1(:,:,i)*Pinv*Winv*A';
            L5 = -A*Winv*Pinv*P1(:,:,i)*Pinv*P1(:,:,j)*Pinv*Winv*A';
            H(i,j) = tIAWinvPinv*(vec(L1+L2+L4+L5));
            H(j,i) = H(i,j);
        end
    end
    H = H + kappa*diag(1./(z.^2) + 1./((1-z).^2)) ;    
    %toc
    %tic
    %H = sigma22.*(2*Pd+V).*V + kappa*diag(1./(z.^2) + 1./((1-z).^2)) ;    
    %[VV,WW]=eig(H);
    %VV
    %diag(WW)
    RR = chol(H);
    Hinvg = (RR\(RR'\g));
    Hinv1 = (RR\(RR'\ones_m));
    dz = -Hinvg + ((ones_m'*Hinvg) / (ones_m'*Hinv1))*Hinv1;
    %toc
    %tic
    deczi = find(dz < 0);
    inczi = find(dz > 0);
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);

    while (1)
        zp = z + s*dz;
        %[P,~,~]=idare(A',U'*diag(sqrt(zp)),Q,R,S,I);
        [~,Ssvd,Vsvd]=svd(diag(sqrt(zp))*U,'econ');
        [P,~,~]=idare(A',Vsvd*Ssvd,Q,sigma2*eye(n),zeros(n,n),I);
        fzp = trace(P) - kappa*sum(log(zp) + log(1-zp));      
        if (fzp <= fz + alpha*s*g'*dz)
            break;
        end
        s = beta*s;
    end
    %toc
    z = zp; fz = fzp; 
    zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
    %P  = diag(zhat./sigma2) - diag(zhat./sigma2)*Ur*(inv(invSr2 + Ur'*diag(zhat./sigma2)*Ur))*Ur'*diag(zhat./sigma2)'; 
    fzhat= fzp;
    zhat=0;
    idata=idata+1;
    time10=toc;
    %[P,~,~]=idare(A',U'*diag(sqrt(zp)),Q,R,S,I);
    [~,Ssvd,Vsvd]=svd(diag(sqrt(zp))*U,'econ');
    [P,~,~]=idare(A',Vsvd*Ssvd,Q,sigma2*eye(n),zeros(n,n),I);
    fz = trace(P) - kappa*sum(log(zp) + log(1-zp));
    fprintf('%4d %10.3f %10.3f %10.3f %10.3f %10.3f \n', iter, s, -g'*dz/2, -fz, trace(P) ,time10);
    data(idata,:)=[iter  -fz trace(P) -fzhat time10];
    if(-g'*dz/2 <= NT_TOL)
        break;
    end
end

[zsort,sensors]=sort(z,'descend');
thres=zsort(m-k); zhat=(z>thres);
%[P,~,~]=idare(A',U'*diag((zhat)),Q,R,S,I);
%[P,~,~]=idare(A',U'*diag((zhat)),Q,R,S,I);
[~,Ssvd,Vsvd]=svd(diag(zhat)*U,'econ');
[P,~,~]=idare(A',Vsvd*Ssvd,Q,sigma2*eye(n),zeros(n,n),I);
L = trace(P);
ztilde = z; 
Utilde = trace(P) + 2*m*kappa;
%z


