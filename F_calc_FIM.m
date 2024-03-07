function [FIM] = F_calc_FIM(U,sensors,Rinv)
    C = U(sensors,:);
    FIM = C'*Rinv*C;
end