function[isensors, det_test]= F_subsensor_DG_p(U, r, p, s, H, isensors, det_test)
%
% load Psi.mat
% load Xorg_sst.mat
% load Xorg_icec.mat
% Xorg=[Xorg_sst;Xorg_icec];
% r=10;
% U=Psi(:,1:r);
%
[n,~]=size(U);
%
% p=8;
% s=2;
% TT=52*10;
% isensors=[9660;4362;7864;7940;4092];
%
% [H]=makesensor_matrix(s, r/s, n, isensors);

[ps,~]=size(isensors);
% ip=zeros(1,p);
%
% ip(1,1:ps)=isensors(1:ps,1)';
ip=isensors';
C=H*U;
CTC=C'*C;

arg=zeros(1,n/s);

for i=1:s
    for pp=1:ps
        U(ip(pp)+(i-1)*n/s,:)=zeros(1,r);
    end
end
%writematrix(U,'U_5.xls')
%csvwrite('U_5.txt',U);
%  save('U_5.csv','U');
C=[];
for i=1:s
    C(:,:,i)=U((1+(i-1)*n/s):(i*n/s) ,1:end);
end
CTCI=inv(CTC);
for pp=(ps+1):p
%    size(CTCI)
     determinant(pp)=det(CTCI);
     [~,e] = eig(CTCI);
     [eigen,~] = sort(diag(e));
%      eigen(1)
%      eigen(2)
%      eigen(3)
%      eigen(4)
%      eigen(5)
%      eigen(6)
    for nn=1:n/s
        arg(nn)=1;
        for ss=1:s
            utmp(ss,:)=C(nn,:,ss);
        end
        
        arg(nn)=det(eye(s,s) +utmp(1:s,:)*CTCI*utmp(1:s,:)');
%         if arg(nn) == 1
%             arg(nn)=0/0;
%         end

%arg(nn)=det(utmp(1:s,:)*CTCI*utmp(1:s,:)');
        
        %        for ss=1:s
        %            abs(nn)=abs(nn)*(1+utmp(ss,:)*CTCI*utmp(ss,:)');
        %            CC=utmp(ss,:)/norm(utmp(ss,:));
        %             for st=ss+1:s
        %                 utmp(st,:)=(utmp(st,:))-(CC*utmp(st,:)')*(CC);
        %             end
        %       end
    end
    [value,ip(pp)]=max(arg);
    det_test(pp)=max(arg);

    for ss=1:s
        % ss
        %ip(pp)+(ss-1)*n/s
        utmp(ss,:)=C(ip(pp),:,ss);
    end
    %     utmp'
    CTCI=CTCI*(eye(r,r) - utmp'*inv(eye(s,s)+utmp*CTCI*utmp')*utmp*CTCI);
    %size(utmp*CTCI*utmp')
    %size(utmp'*inv(eye(s,s)+utmp*CTCI*utmp')*utmp*CTCI)

    for st=1:s
        C(ip(pp),:,st)=zeros(1,r);
    end
    
    %                for ss=1:s
    %                 CC=C(ip(pp),:,ss);
    %                 CC=CC/norm(CC);
    %                 for st=1:s
    %                     C(:,:,st)=C(:,:,st)-C(:,:,st)*CC'*CC;
    %                 end
    %             end
    
    
    %    save('U_',(pp),'.csv','C');
    %    csvwrite('U_',pp,'.txt', C);
    %    C(n+ip(pp),:)=zeros(1,r);
    CTC=CTC+utmp'*utmp;
    inv(CTC);
    isensors=ip';
    %size(CTC)
    %     [H]=makesensor_matrix(s, pp, n, isensors);
    %     [Xestimate_DG]= reconst(Xorg,H,Psi,r);
    %     [Perror_DG,Perror_std_DG] = subcalculation_error(TT,Xorg, Psi, Xestimate_DG,r);
    %     Perror_DG
end
%save('determinant.mat','determinant');
end

%     for ss=1:s
%         utmp(ss,:)=U(n/s*(ss-1)+ip(pp),:);
%     end
%     for ss=1:s;
%         U(n/s*(ss-1)+ip(pp),:)=zeros(1,r);
%     end
%     det(CTC);
%     det(inv(CTCI));


% for i=1:n
%     rownorm(i)=norm(U(i ,:));
% end
% rownorm'
%  save('rownorm_U.mat','rownorm');
% aaa