function[isensors, det_test]= F_subsensor_DG_r(U, r, p, s, det_test)
% s=2;
% p=5;
% load Psi_PIV.mat;
% load S_PIV.mat;
% load V_PIV.mat;
% Psi=Psi_PIV;
% S=S_PIV;
% V=V_PIV;
% Xorg=Psi*S*V';
%
% Xorg_u=Xorg(1:3048,:);
% Xorg_v=Xorg(3049:6096,:);
%
% [Psi_u, S_u, V_u] = svd(Xorg_u,'econ');
% [Psi_v, S_v, V_v] = svd(Xorg_v,'econ');
%
% size(Psi_u);
% size(Psi_v);
%
% %U=[Psi_u(:,1:50) Psi_v(:,1:5)];
% r1=50;
% r2=5;
%
% %U=Psi_u(:,1:r1);
% %V=Psi_v(:,1:r2);
%
% [n,~]=size(Psi_u);
% U1=[Psi_u(:,1:r1), zeros(n,r2)];
% U2=[zeros(n,r1),Psi_v(:,1:r2)];
%
% U=[U1;U2];
[n,~]=size(U);

C=[];
for i=1:s
    C(:,:,i)=U((1+(i-1)*n/s):(i*n/s) ,1:end);
end

THETApp=zeros(p*s,r);
THETAppTHETApptinv=zeros(p*s,p*s);
utmp=zeros(s,r);
abs=zeros(1,n/s);
sensor=zeros(1,p);

for pp=1:p
    Y=eye(r,r)-THETApp(1:(pp-1)*s,:)'*THETAppTHETApptinv(1:(pp-1)*s,1:(pp-1)*s)*THETApp(1:(pp-1)*s,:);
    for nn=1:n/s
        abs(nn)=1;
        for ss=1:s
            utmp(ss,:)=C(nn,:,ss);
        end
        abs(nn)=det(utmp(1:s,:)*Y*utmp(1:s,:)');
        
    end
    %    save('abs_u.mat','abs');
    [~,sensor(pp)]=max(abs);
    det_test(pp)=max(abs);
    
    for ss=1:s
        THETApp(pp*s-(s-ss),:)=C(sensor(pp),:,ss);
        THETAppTHETApptinv(1:pp*s-(s-ss),1:pp*s-(s-ss))=inv(THETApp(1:pp*s-(s-ss),:)*THETApp(1:pp*s-(s-ss),:)');
    end
    
    for st=1:s
        C(sensor(pp),:,st)=zeros(1,r);
    end
    
%         for ss=1:s
%             CC=C(sensor(pp),:,ss);
%             CC=CC/norm(CC);
%             for st=1:s
%                 C(:,:,st)=C(:,:,st)-C(:,:,st)*CC'*CC;
%             end
%         end
end
isensors=sensor';
end

%         for ss=1:s
%             abs(nn)=abs(nn)*utmp(ss,:)*Y*utmp(ss,:)';
%             CC=utmp(ss,:)/norm(utmp(ss,:));
%             for st=ss+1:s
%                 utmp(st,:)=(utmp(st,:))-(CC*utmp(st,:)')*(CC);
%             end
%         end


%         for ss=1:s
%         CC=C(sensor(pp),:,ss);
%         CC=CC/norm(CC);
%         for st=1:s
%             C(:,:,st)=C(:,:,st)-C(:,:,st)*CC'*CC;
%         end
%     end