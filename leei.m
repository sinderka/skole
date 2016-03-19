function [utdata] = leei()
m = 50;
k = 50;
n = 1; %%%%%%må endres på!!
func =2; % må være 3
utdata = zeros(1,3);

Y = linspace(0,1,m+2)'; X = Y;
%hh = Y(2)-Y(1);
hs = 1/(m+1);
A = -1/hs^2*gallery('poisson', m);
%A = 1*hs^2*gallery('tridiag', m^2);
%A = hs^4*eye(m^2);
%A = zeros(m^2);
T = linspace(0,1,k); ht = T(2)-T(1);
intmet = 2;


if func == 1
    solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1); 
    f1 = @(t) 1/(t+1).^2;   g1 = @(x,y) x*(x-1)*y.*(y-1);
    f2 = @(t) -t/(t+1);     g2 = @(x,y) 2*x.*(x-1)+2*y.*(y-1);
elseif func == 2
    solution1 = @(t,x,y) sin(x*y).*(x-1).*(y-1).*t;
    f1 = @(t) 1;            g1 = @(x,y) sin(x.*y).*(x-1).*(y-1);
    f2 = @(t) -t;           g2 = @(x,y) (y-1).*y.*(2*cos(x.*y)-(x-1).*y.*sin(x*y)) + (x-1).*x.*(2*cos(x.*y)-(y-1).*x.*sin(x*y));
elseif func == 3
    solution1 = @(t,x,y)sin(pi*x).*sin(pi*y).*sin(t);
    g1 = @(x,y) sin(pi*x).*sin(pi*y);f1 = @(t) cos(t)+2*pi^2*sin(t);
    f2 = @(t) t-t ; g2 = @(x,y) x;
end

sol1 = zeros(m^2,k);
for i = 1:k
    for j = 1:m
        sol1(1+(j-1)*(m):j*(m),i) = solution1(T(i),X(j+1),Y(2:end-1));
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
[U,ant1] = locprojmeth(F1,A,v1,k,ht,n,intmet);
[Utemp,ant2] = locprojmeth(F2,A,v2,k,ht,n,intmet);
utdata(1) = max(ant1,ant2);
U = Utemp + U;
U1 = locprojmeth1(F1,A,v1,k,ht,m^2,intmet);
U1 = U1+ locprojmeth1(F2,A,v2,k,ht,m^2,intmet);
%U = zeros(m^2,k);
%U1 = U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% alternativ måte %%%
%%% intigrer: du/dt = 1/h^2*A*u(t)+f(t)
%ZZZ = zeros(m^2);

%U2 = v1*F1;


%U2 = v1*F1-U2;

U2 = integrere(m,k,v1*F1,ht,A);
U2 = U2 + integrere(m,k,v2*F2,ht,A);
%Utest = A\U2;
%U2  = A*U2;
%U2 = v1*F1-U2;
%figure(100) 
%plot(Utest)




%%% alternativ måte %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
utdata(1,2) = toc;
Err = zeros(1,m^2);
for i = 1:m^2
    Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
end

utdata(1,3) = max(Err);

figure(1)
plot(U); title('krylov, restart')
figure(2)
plot(sol1); title('rett')
figure(3)
plot(sol1-U); title('rett-krylov, restart')
figure(4)
plot(U1); title('krylov, ikke restart')
figure(5)
plot(sol1-U1); title('rett-krylov, ikke restart')
figure(6)
plot(U2); title('direkte')
figure(7)
plot(sol1-U2); title('rett-direkte')

end

function [U,ant] = locprojmeth(Zn,A,v,k,ht,n,intmet)
U = zeros(length(A),k);
diff =10;
hn = norm(v,2);
ant = 0;
while diff > 10^-16
    [V,H,hm] = krylov(A,v,n);
    [Zn] = locintmeth(Zn(end,:),H,hn,k,ht,n,intmet);
    hn = hm;
    ns = V(:,1:n)*Zn;
    U = U + ns;
    v = V(:,end);
    diff = max(max(abs(ns)));
    ant = ant+1;
end
end

function [U] = locprojmeth1(Zn,A,v,k,ht,n,intmet)
hn = norm(v,2);
[V,H,~] = krylov(A,v,n);
[Zn] = locintmeth(Zn(end,:),H,hn,k,ht,n,intmet);
U = V(:,1:n)*Zn;
end



function [Zn] = locintmeth(F,H,hn,k,ht,n,intmet)
%%%%%Solves the equation z'-H*z=hn*F Numerically
%%% ht is stepsize, m is the size of the square matrix H
%%% len_T is the number of time-steps, zn is a vector of length len_T
%%% intmet says something about which integration-method that should be
%%% used. Boundary equal to zero.

if intmet == 1
    Zn = zeros(n,k);
    for j = 2:1:k
        Zn(:,j) = Zn(:,j-1) + H*ht*Zn(:,j-1);
        Zn(1,j) = Zn(1,j) + ht*hn*F(j-1);
    end
elseif intmet == 2
    mat = inv(eye(n)-ht/2*H);
    e1 = zeros(n,1); e1(1) = 1;
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*(F(1)*hn+F(2)*hn))*e1;
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*H*Zn(:,j-1))+mat*(ht/2*(F(j-1)*hn+F(j)*hn))*e1;
        %Zn(1,j) = Zn(1,j) 
    end
else
    disp('error in "intmeth"');
end
end

function U = integrere(m,k,F,ht,A)
U = zeros(m^2,k);
mat = inv(eye(m^2)-(ht/2)*A);
U(:,2) = mat*(ht/2*(F(:,2)+F(:,1)));
%U(:,2) = (ht/2*(F(:,2)+F(:,1)));
for j = 3:k
    U(:,j) = mat*(U(:,j-1) +(ht/2*(F(:,j)+F(:,j-1)+A*U(:,j-1))));
    %U(:,j) = (U(:,j-1) +(ht/2*(F(:,j)+F(:,j-1)+A*U(:,j-1))));
end
end