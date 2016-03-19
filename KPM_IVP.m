function [z] = KPM_IVP(H,T,len_T,Hn,zn,M)
    z = zeros(M,len_T);
        for j = 2:1:len_T
             h = T(j)-T(j-1);
             z(:,j) = z(:,j-1) + h*(H*z(:,j-1));
             z(1,j) = z(1,j) + h*Hn*zn(j-1);
        end
end