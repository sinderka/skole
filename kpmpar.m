function [ data,err ] = kpmpar( k,m,n,prob,para,numf )
% k = number of time steps
% m = number of spacial steps in each direction
% n = restartvariabel
% prob = 1: solve problem 1
%        2: solve problem 2
%%%%%% Ordner returdata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = zeros(1,5);      % data skal inneholde k n m max(Err) tid
data(1,1) = k; data(1,2) = m; data(1,3) = n; err =0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic;
parpool(para);
hs = 1/(m+1);
A = hs^2*gallery('poisson', m);
T = linspace(0,1,k); ht = T(2)-T(1);
intmet = 2;
Y = linspace(0,1,m+2)'; X = Y;
U=0;

if prob == 1
    
    solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1)*numf;
    sol1 = zeros(m^2,k);
    for i = 1:k
        for j = 1:m
            sol1(1+(j-1)*(m):j*(m),i) = solution1(T(i),X(j+1),Y(2:end-1));
        end
    end
    
    f1 = @(t) 1/(t+1).^2;
    g1 = @(x,y) x*(x-1)*y.*(y-1);
    
    v1 = zeros(m^2,1);
    for i = 1:m
        v1(1+(i-1)*m:i*m) = g1(X(i+1),Y(2:m+1));
    end
    F1 = zeros(1,k);
    for i = 1:k
        F1(i) = f1(T(i));
    end
    tic;
    parfor i = 1:numf
    U = U + projmeth(F1,A,v1,k,ht,n,intmet);
    end
    data(1,5) = toc;
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
    end
    
    data(1,4) = max(Err);
end
if prob == 2
    
    solution2 = @ (t,x,y) t*x*(x-1)*sin(pi*y)*numf;
    sol2 = zeros(m^2,k);
    for i = 1:k
        for j = 1:m
            sol2(1+(j-1)*(m):j*(m),i) = solution2(T(i),X(j+1),Y(2:end-1));
        end
    end
    
    f2 = @(t) t-t+1;
    g2 = @(x,y) x*(x-1)*sin(pi*y);
    v2 = zeros(m^2,1);
    for i = 1:m
        v2(1+(i-1)*m:i*m) = g2(X(i+1),Y(2:m+1));
    end
    F2 = zeros(1,k);
    for i = 1:k;
        F2(i) = f2(T(i))';
    end
    tic;
    parfor i = 1:numf
    U = U+ projmeth(F2,A,v2,k,ht,n,intmet);
    end
    data(1,5) = toc;
    
    
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol2(i,:)-U(i,:) ,Inf)/max(abs(sol2(i,:)));
    end
    
    data(1,4) = max(Err);
end

delete(gcp);
end

function [U,err] = projmethtest(Zn,A,v,k,ht,m,intmet,sol)
err = zeros(1,30);
i = 0;
U = zeros(length(A),k);
diff =1;
hn = norm(v,2);
while diff > 10^(-15) && i <= 30
    [V,H,hm] = krylov(A,v,m);
    [Zn] = intmeth(Zn(end,:),H,hn,k,ht,m,intmet);
    hn = hm;
    ns = V(:,1:m)*Zn;
    U = U + ns;
    v = V(:,end);
    diff = max(max(abs(ns)));
    i = i+1;
    err(i) = max(max(abs(sol-U)));
end

end



