function [U] = projmeth(Zn,A,v,k,ht,n,intmet)
U = zeros(length(A),k);
diff =1;
hn = norm(v,2);
while diff > 10^(-15)
    [V,H,hm] = krylov(A,v,n);
    [Zn] = intmeth(Zn(end,:),H,hn,k,ht,n,intmet);
    hn = hm;
    ns = V(:,1:n)*Zn;
    U = U + ns;
    v = V(:,end);
    diff = max(max(abs(ns)));
end