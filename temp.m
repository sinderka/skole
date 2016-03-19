function utdata = temp
m = 3; n = m^2; k=10;
utdata = zeros(1,3);
%err = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic;
hs = 1/(m+1);
A = 1/hs^2*gallery('poisson', m);
%A = hs^2*eye(m^2);
%A = hs^2*gallery('poisson', m);
T = linspace(0,1,k); ht = T(2)-T(1);
intmet = 2;
Y = linspace(0,1,m+2)'; X = Y;
U = 0;
%    solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1); f1 = @(t) 1/(t+1).^2; g1 = @(x,y) x*(x-1)*y.*(y-1);
%    solution1 = @(t,x,y) exp(x.*y).*sin(pi*x).*y.*(y-1).*cos(t)*t; f1 = @(t) cos(t)-t.*sin(t); g1 = @(x,y) exp(x.*y).*sin(pi*x).*y.*(y-1);
solution = @(t,x,y) sin(pi*x)*sin(pi*y)*sin(t);
%g1 = @(x,y) sin(pi*x)*sin(pi*y); f1 = @(t) sin(t);
%g2 = @(x,y) -2*pi^2*sin(pi*x)*sin(pi*y); f2 = @(t) cos(t);
g1 = @(x,y) x.*(x-1).*y.*(y-1); f1 = @(t) 1;
g2 = @(x,y) 2*y.*(y-1) +x.*(x-1); f2 = @(t) 1;

sol = zeros(m^2,k);
for i = 1:k
    for j = 1:m
        sol(1+(j-1)*(m):j*(m),i) = solution(T(i),X(j+1),Y(2:end-1));
    end
end


v1 = zeros(m^2,1);
v2 = zeros(m^2,1);
for i = 1:m
    v1(1+(i-1)*m:i*m) = g1(X(i+1),Y(2:m+1));
    v2(1+(i-1)*m:i*m) = g2(X(i+1),Y(2:m+1));
end
F1 = zeros(1,k);
F2 = zeros(1,k);
for i = 1:k
    F1(i) = f1(T(i));
    F2(i) = f2(T(i));
end

tic;
U = U + -locprojmeth(F1,A,v2,k,ht,n,intmet);
%U = U + -locprojmeth(F2,A,v2,k,ht,n,intmet);
utdata(1,2) = toc;
Err = zeros(1,m^2);
%for i = 1:m^2
%    Err(i) = norm(sol(i,:)-U(i,:) ,Inf)/max(abs(sol(i,:)));
%end

utdata(1,3) = max(Err);
%figure(1)
%plot(U)
%figure(2)
%plot(sol)
%figure(3)
%plot(sol-U)
%figure(4)
%plot(sol/max(max(abs(sol)))-U/max(max(abs(U))))

x = A\v2;
figure(4)
plot(v1); title('rett løsning')


figure(5)
plot(-x); title('direkte løsning')

figure(6)
plot(U); title('krylov')
figure(7)
plot(x+U); title('rett - krylov')
figure(8)
plot(v1+x); title('rett - direkte')
figure(9)
plot(v1-U); title('direkte - krylov')


end


function [U] = locprojmeth(Zn,A,v,k,ht,n,intmet)
%U = zeros(length(A),1);
%diff =1;
hn = norm(v,2);
%while diff > 10^(-5)
    [V,H,hm] = krylov(A,v,n);
    %[V,H] = Arnoldi(A,v,n);
    %hm = H(n+1,n);
    H = H(1:n,1:n);
    U = hn*V(:,1:n)*(H\(V(:,1:n)'*V(:,1)));
    A-V(:,1:n)'*(H*(V(:,1:n)))
    %hn = hm;
    %U = U +ns;
    %v = V(:,end);
    %diff = max(max(abs(ns)));
%end
end