function [zhat L ztilde Utilde,sensors, data] = F_sens_sel_A_approxnt_corr(A, Ur, Sr2, sigma2, k, MAXITER)
% Solves the problem
%	maximize log det (sum_{i=1}^m z_i a_i a_i^T) + kappa sum_{i=1}^m(log(z_i)+ log(1-z_i))
%	subject to sum(z) = k
%			   0 <= z_i <= 1, i=1,..., m
% variable z in R^m
% problem parameters kappa (>0), a_1, ..., a_m in R^n
%
% see paper Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% Nov 2007 Siddharth Joshi & Stephen Boyd
sigma2=sigma2;
% Newton's method parameters
%MAXITER  = 30;
NT_TOL = 1e-7;
GAP = 1.005;
GAP = 1.001;
%GAP = 1.0001;
% Backtracking line search parameters
alpha = 0.01;
beta = 0.5;

[m n] = size(A);
z = ones(m,1)*(k/m); % initialize
g = zeros(m,1);
ones_m = ones(m,1);
kappa = log(GAP)*n/m; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

fprintf('\nIter.  Step_size  Newton_decr.  Objective  trinv trinv_dis time \n');
invSr2=inv(Sr2);
P  = diag(z./sigma2) - diag(z./sigma2)*Ur*(inv(invSr2 + Ur'*diag(z./sigma2)*Ur))*Ur'*diag(z./sigma2); 
fz = trace(inv(A'*P*A)) - kappa*sum(log(z) + log(1-z));

fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, trace(inv(A'*P*A)));
tic
idata=0;
for iter=1:MAXITER
    P  = inv(Ur*Sr2*Ur'+diag(sigma2./z));
    %P  = diag(z./sigma2) - diag(z./sigma2)*Ur*(inv(invSr2 + Ur'*diag(z./sigma2)*Ur))*Ur'*diag(z./sigma2); 
    %Pd =                   diag(z./sigma2)*Ur*(inv(invSr2 - Ur'*diag(z./sigma2)*Ur))*Ur'*diag(z./sigma2)'; 
    W = inv(A'*P*A);
    V = P*A*W*W*A'*P;
    Vs = P*A*W*A'*P;

    g = -diag(sigma2./(z.^2))*diag(V)- kappa*(1./z - 1./(1-z));
    sigma22=((sigma2./z.^2)*(sigma2./z.^2)');
    H = -2*sigma22.*P.*V + 2*sigma22.*Vs.*V + 2*diag(sigma2./(z.^3)).*V + kappa*diag(1./(z.^2) + 1./((1-z).^2)) ;    
    %H = sigma22.*(2*Pd+V).*V + kappa*diag(1./(z.^2) + 1./((1-z).^2)) ;    
    %[VV,WW]=eig(H);
    %VV
    %diag(WW)
    R = chol(H);
    Hinvg = (R\(R'\g));
    Hinv1 = (R\(R'\ones_m));
    dz = -Hinvg + ((ones_m'*Hinvg) / (ones_m'*Hinv1))*Hinv1;

    deczi = find(dz < 0);
    inczi = find(dz > 0);
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);

    while (1)
        zp = z + s*dz;
        P  = inv(Ur*Sr2*Ur'+diag(sigma2./zp));
        %P  = diag(zp./sigma2) - diag(zp./sigma2)*Ur*(inv(invSr2 + Ur'*diag(zp./sigma2)*Ur))*Ur'*diag(zp./sigma2)';                    
        fzp = trace(inv(A'*P*A)) - kappa*sum(log(zp) + log(1-zp));

        if (fzp <= fz + alpha*s*g'*dz)
            break;
        end
        s = beta*s;
    end
    z = zp; fz = fzp; 
    zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
    %P  = diag(zhat./sigma2) - diag(zhat./sigma2)*Ur*(inv(invSr2 + Ur'*diag(zhat./sigma2)*Ur))*Ur'*diag(zhat./sigma2)'; 
    fzhat= fzp;
    zhat=0;
    idata=idata+1;
    time10=toc;
    P  = inv(Ur*Sr2*Ur'+diag(sigma2./zp));
    %P  = diag(z./sigma2) - diag(z./sigma2)*Ur*(inv(invSr2 + Ur'*diag(z./sigma2)*Ur))*Ur'*diag(z./sigma2)';    
    APA=(A'*P*A);
    fprintf('%4d %10.3f %10.3f %10.3f %10.3f %10.3f \n', iter, s, -g'*dz/2, -fz, trace(inv(APA)) ,time10);
    data(idata,:)=[iter  -fz trace(inv(APA)) -fzhat time10];
    if(-g'*dz/2 <= NT_TOL)
        break;
    end
end

[zsort,sensors]=sort(z,'descend');
thres=zsort(m-k); zhat=(z>thres);
P  = inv(Ur*Sr2*Ur'+diag(sigma2./zhat));
%P  = diag(zhat./sigma2) - diag(zhat./sigma2)*Ur*(inv(invSr2 + Ur'*diag(zhat./sigma2)*Ur))*Ur'*diag(zhat./sigma2)'; 
L = trace(inv(A'*P*A));
ztilde = z; 
Utilde = trace(inv(A'*P*A)) + 2*m*kappa;
z


