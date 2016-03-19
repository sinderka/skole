function [ inndata,utdata ] = kpmtester( k,m,n,prob,para,func)
% k = number of time steps
% m = number of spacial steps in each direction
% n = restartvariabel
% prob = 1: solve problem 1
%        2: solve problem 1, with a projection method that compares the
%        result to the solution each iteration
%%%%%% Ordner returdata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inndata = zeros(1,6);      % data skal inneholde k n m max(Err) tid
inndata(1,1) = k; inndata(1,2) = m; inndata(3) = n; inndata(1,4) = prob;  inndata(1,5) = para; inndata(1,6) = func;
utdata = zeros(1,3);
%err = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic;
hs = 1/(m+1);
A = hs^2*gallery('poisson', m);
T = linspace(0,1,k); ht = T(2)-T(1);
intmet = 2;
Y = linspace(0,1,m+2)'; X = Y;
U=0;

if func == 1
    solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1); f1 = @(t) 1/(t+1).^2; g1 = @(x,y) x*(x-1)*y.*(y-1);
elseif func == 2
    solution1 = @(t,x,y) exp(x.*y).*sin(pi*x).*y.*(y-1).*cos(t)*t; f1 = @(t) cos(t)-t.*sin(t); g1 = @(x,y) exp(x.*y).*sin(pi*x).*y.*(y-1);
end

sol1 = zeros(m^2,k);
for i = 1:k
    for j = 1:m
        sol1(1+(j-1)*(m):j*(m),i) = solution1(T(i),X(j+1),Y(2:end-1));
    end
end


v1 = zeros(m^2,1);
for i = 1:m
    v1(1+(i-1)*m:i*m) = g1(X(i+1),Y(2:m+1));
end
F1 = zeros(1,k);
for i = 1:k
    F1(i) = f1(T(i));
end

if prob == 1        % projektion method
    tic;
    [U,utdata(1)] = projmeth1(F1,A,v1,k,ht,n,intmet);
    utdata(1,2) = toc;
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
    end
    
    utdata(1,3) = max(Err);
    inndata(1,5) = 0;
end

if prob == 2        % parallell projektion method
    ant = 0;
    parpool(para);
    tic;
    parfor i = 1:m^2
        e = zeros(m^2,1); e(i) = v1(i);
        [U1,ant1] = projmeth1(F1,A,e,k,ht,n,intmet);
        U = U + U1;
        ant = max(ant,ant1);
    end
    utdata(1,2) = toc;
    delete(gcp)
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
    end
    utdata(1) = ant;
    utdata(1,3) = max(Err);
    %data(1,3) = n;
end

if prob == 3        % ikke projektion method
    tic;
    U = intmeth2(F1,A,v1,k,ht,m^2,intmet);
    utdata(1,2) = toc;
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
    end
    utdata(1,3) = max(Err);
    utdata(1,1) = 1;
    inndata(1,5) = 0;
end

if prob == 4        % parallell ikke projektion method
    parpool(para);
    tic;
    parfor i = 1:m^2
        e = zeros(m^2,1); e(i) = v1(i);
        U = U + intmeth2(F1,A,e,k,ht,m^2,intmet);
    end
    utdata(1,2) = toc;
    delete(gcp)
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
    end
    utdata(1,3) = max(Err);
    utdata(1) = 1;
end


if prob == 5 %(ikke bruk uten at du vet hva du gjør)
    tic;
    [U1,err] = projmethtest(F1,A,v1,k,ht,n,intmet,sol1);
    U = U + U1;
    data(1,5) = toc;
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
    end
    
    data(1,4) = max(Err);
    data(1,3) = n;
    data(1,8) = 0;
end

end

function [U,ant] = projmethtest(Zn,A,v,k,ht,m,intmet)
%err = zeros(1,30);
%i = 0;
U = zeros(length(A),k);
diff =1;
hn = norm(v,2);
ant = 0;
while diff > 10^(-15) %&& i <= 30
    [V,H,hm] = krylov(A,v,m);
    [Zn] = intmeth(Zn(end,:),H,hn,k,ht,m,intmet);
    hn = hm;
    ns = V(:,1:m)*Zn;
    U = U + ns;
    v = V(:,end);
    diff = max(max(abs(ns)));
    ant = ant+1;
    %err(i) = max(max(abs(sol-U)));
end

end
function [Zn] = intmeth2(F,A,v,k,ht,m,intmet)
%%%%%Solves the equation z'(t)-A*z(t)=F(t)*v Numerically
%%% ht is stepsize, m is the size of the square matrix H
%%% len_T is the number of time-steps, zn is a vector of length len_T
%%% intmet says something about which integration-method that should be
%%% used. Boundary equal to zero.

if intmet == 1
    Zn = zeros(m,k);
    for j = 2:1:k
        Zn(:,j) = Zn(:,j-1) + A*ht*Zn(:,j-1)+ht*F(j-1)*v;
    end
elseif intmet == 2
    Zn = zeros(m,k);
    Zn(:,2) = ht/2*(F(1)*v+F(2)*v);
    for j = 3:1:k
        Zn(:,j) = Zn(:,j-1)+ht/2*A*Zn(:,j-1) +ht/2*(F(j-1)*v+F(j)*v) ;
    end
else
    disp('error in "intmeth"');
end
end

function [U,ant] = projmeth1(Zn,A,v,k,ht,n,intmet)
U = zeros(length(A),k);
diff =1;
hn = norm(v,2);
ant = 0;
while diff > 10^(-15)
    [V,H,hm] = krylov(A,v,n);
    [Zn] = intmeth(Zn(end,:),H,hn,k,ht,n,intmet);
    hn = hm;
    ns = V(:,1:n)*Zn;
    U = U + ns;
    v = V(:,end);
    diff = max(max(abs(ns)));
    ant = ant+1;
end
end





