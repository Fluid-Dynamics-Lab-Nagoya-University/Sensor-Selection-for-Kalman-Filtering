% function to select sensor based on Gramian
% F_sensor_Gram.m : ----------- Programer : Kumi Nakai, Taku Nonomura 2020
%                               Last modified: 2021/11/17 K.Nakai(maintanence)

function [time, H, isensor]=F_sensor_Gram(Aorg,Corg,p)

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
            Dprev=[];
            H=zeros(1,n);
        end
        maxrank=0;
        for i=1:1:n
            C=[Cprev;Ctmp(i,:)];
            CTC = C'*C;
            vecCTC = reshape(CTC,[],1);
            vecWo = - AAIinv * vecCTC; % fastest calculation (200726)
            G(:,:,i) = reshape(vecWo,[r,r]);
            ranki(i) = rank(G(:,:,i));
            maxrank = max(rank(G(:,:,i)),maxrank);
        end
        for i=1:1:n
            if( ranki(i)==maxrank )
                Gram=G(:,:,i);
                obj1(i) = maxrank;
                obj2(i) = det(Gram);
            end
        end

        [~,isensor(pp)] = max(obj2);
        H(pp,isensor(pp)) = 1;
        text = [ 'Result(Gramian-based): isensor=' num2str(isensor(pp)) ', rank(Wo)=' num2str(obj1(isensor(pp))) ];
        disp ( text ) ;
        Ctmp(isensor(pp),:) = zeros(1,r); %zero reset
        
    end
    time=toc;
    isensor=isensor';
end

