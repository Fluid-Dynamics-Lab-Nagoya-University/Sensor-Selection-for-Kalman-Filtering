function [Wo]=F_calc_Gram_vi(A,U,sensors,iR)
[ps,~]=size(sensors);
[n,r]=size(U);
AAI1 = kron(A',A')-eye(r^2);
iAAI1 = inv(AAI1);
C = U(sensors,:);
if isempty(iR)
    iR = eye(ps);
end
vecCC = reshape(C'*iR*C,[],1);  %<------- CdC1
vecWo = - iAAI1 * vecCC;
Wo = reshape(vecWo,[r,r]);
end