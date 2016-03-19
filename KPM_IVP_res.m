function [y] = KPM_IVP_res(A,v,M,zn,h,len_T)
y = 0;
diff =100;
Hn = norm(v,2);
while diff > 0.0001 %%%Finne en bedre måte å rengne ut feilen på?
    [V,H] = Arnoldi(A,v,M);
    %[V,H,Hm] = krylov(A,v,M);
    [zn] = intmeth(zn(end,:),H,Hn,len_T,h,M);
    if (M+1 < length(v))
    Hn = H(M+1,M);  %% For Arnoldi
    else
        Hn = 0;
    end
    %Hn = Hm;        %%  For krylov
    ns = V(:,1:M)*zn;
    y = y + ns;
    v = V(:,end);
    diff = max(max(ns));
end
end