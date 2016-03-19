function [ utdata ] = seriell(m,k,n,prob,func,conv )

%%% Legger til verdier om de ikke er fastsatt
if nargin == 0
    m = 5;
    k = 5;
    n = 1;
    prob = 1;
    func  = 1;
end
%%%

%%% Legger til conv om denne ikke er fastsatt
if nargin < 6
   conv = 10^-12; 
end
%%%

%%% Initsiell data
utdata = zeros(1,3);
Y = linspace(0,1,m+2)'; X = Y;
hs = 1/(m+1);
A = -1/hs^2*gallery('poisson', m);
T = linspace(0,1,k); ht = T(2)-T(1);
U = 0;
ant = 1;
%%%


%%% Velg en funksjon å teste for
if func == 1
    solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1);
    f1 = @(t) 1/(t+1).^2;   g1 = @(x,y) x*(x-1)*y.*(y-1);
    f2 = @(t) -t/(t+1);     g2 = @(x,y) 2*x.*(x-1)+2*y.*(y-1);
elseif func == 2
    solution1 = @(t,x,y) exp(x.*y).*y.*(y-1).*sin(pi*x).*t.*cos(t);
    f1 = @(t) cos(t)-t*sin(t);            g1 = @(x,y) (y-1).*y.*exp(x.*y).*sin(pi*x);
    f2 = @(t) -t.*cos(t);           g2 = @(x,y) exp(x.*y).*(x.^2.*(y-1).*y+x.*(4*y-2)+2).*sin(pi*x) +(y-1).*y.^3.*exp(x.*y).*sin(pi*x)+...
    2*pi*(y-1).*y.^2.*exp(x.*y).*cos(pi*x)-pi^2*(y-1).*y.*exp(x.*y).*sin(pi*x);
end
%%%

%%% Finner løsning, v og F
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
%%%


%%% Løser det faktiske problemet med 
if prob == 1 %% Restarted krylov
    tic;
    [U1,temp1] = locprojmeth(F1,A,v1,k,ht,n,conv);
    [U2,temp2] = locprojmeth(F2,A,v2,k,ht,n,conv);
    U = U1+U2;
    utdata(1) = max(temp1,temp2);
    utdata(2) = toc;
elseif prob == 2 %% Full krylov
    tic;
    U1 = locprojmeth1(F1,A,v1,k,ht,m^2);
    U2 = locprojmeth1(F2,A,v2,k,ht,m^2);
    U = U1+U2;
    utdata(1) = 1;
    utdata(2) = toc;
elseif prob == 3 %% Direkte integrasjon
    tic;
    U = integrere(m,k,v1*F1+v2*F2,ht,A);
    utdata(1) = 1;
    utdata(2) = toc;
end
%%%


%%% Beregner den største relative feilen
Err = zeros(1,m^2);
for i = 1:m^2
    Err(i) = norm(sol1(i,:)-U(i,:) ,Inf)/max(abs(sol1(i,:)));
end
utdata(3) = max(Err);
%%%

end

function U = integrere(m,k,F,ht,A) %%% Intigrere direkte
%%%%%Solves the equation z'-A*z=v*F
%%% ht is stepsize, m^2 is the size of the square matrix A
%%% k is the number of time-steps.
%%% Zn(:,1) = 0.
U = zeros(m^2,k);
mat = inv(eye(m^2)-(ht/2)*A);
U(:,2) = mat*(ht/2*(F(:,2)+F(:,1)));
for j = 3:k
    U(:,j) = mat*(U(:,j-1) +(ht/2*(F(:,j)+F(:,j-1)+A*U(:,j-1))));
end
end
function [Zn] = locintmeth(F,H,hn,k,ht,n) %%% intigrere KPM %%% dette må gjøers på en smartere måte!!!!!!! 
%%%%%Solves the equation z'-H*z=hn*F Numerically
%%% ht is stepsize, n is the size of the square matrix H
%%% k is the number of time-steps.
%%% Zn(:,1) = 0.
    mat = inv(eye(n)-ht/2*H);
    e1 = zeros(n,1); e1(1) = 1;
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*hn*e1*(F(1)+F(2)));
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+hn*e1*(F(j)+F(j-1))));
    end
end
function [U] = locprojmeth1(Zn,A,v,k,ht,n) %%% Full KPM
%%%%%Projection method for the full krylov space
hn = norm(v,2);
[V,H,~] = krylov(A,v,n);
[Zn] = locintmeth(Zn(end,:),H,hn,k,ht,n);
U = V(:,1:n)*Zn;
end
function [U,ant] = locprojmeth(Zn,A,v,k,ht,n,conv) %%% Restarted KPM
%%%%%Projection method for any kyylov space
U = zeros(length(A),k);
diff =10;
hn = norm(v,2);
ant = 0;
while diff > conv
    [V,H,hm] = krylov(A,v,n);
    [Zn] = locintmeth(Zn(end,:),H,hn,k,ht,n);
    hn = hm;
    ns = V(:,1:n)*Zn;
    U = U + ns;
    v = V(:,end);
    diff = max(max(abs(ns)));
    ant = ant+1;
end
end

















