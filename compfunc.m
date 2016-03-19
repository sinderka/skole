function [ data ] = compfunc( k,m,prob )

%%%%%% Ordner returdata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = zeros(1,5);      % data skal inneholde k n m max(Err) tid
data(1,1) = k; data(1,2) = m; data(1,3) = m^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hs = 1/(m+1);
A = hs^2*gallery('poisson', m);
T = linspace(0,1,k); ht = T(2)-T(1);
intmet = 2;
Y = linspace(0,1,m+2)'; X = Y;
U=0;




if prob == 1
    f1 = @(t) 1/(t+1).^2;
    
    Y = linspace(0,1,m+2)'; X = Y;
    
    g1 = @(x,y) x*(x-1)*y.*(y-1);
    
    
    v1 = zeros(m^2,1);
    for i = 1:m
        v1(1+(i-1)*m:i*m) = g1(X(i+1),Y(2:m+1));
    end
    
    F1 = zeros(1,k);
    for i = 1:k
        F1(i) = f1(T(i));
    end
    
    tic; U = U + intmeth2(F1,A,v1,k,ht,m^2,intmet); data(1,5) = toc;
    
    solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1);
    sol2 = zeros(m^2,k);
    for i = 1:k
        for j = 1:m
            sol2(1+(j-1)*(m):j*(m),i) = solution1(T(i),X(j+1),Y(2:end-1));
        end
    end
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol2(i,:)-U(i,:) ,Inf)/max(abs(sol2(i,:)));
    end
    data(1,4) = max(Err);
end

if prob == 2
    f2 = @(t) t-t+1;
    
    Y = linspace(0,1,m+2)'; X = Y;
    
    g2 = @(x,y) x*(x-1)*sin(pi*y);
    
    
    v2 = zeros(m^2,1);
    for i = 1:m
        v2(1+(i-1)*m:i*m) = g2(X(i+1),Y(2:m+1));
    end
    
    F2 = zeros(1,k);
    for i = 1:k
        F2(i) = f2(T(i));
    end
    
    tic; U = U + intmeth2(F2,A,v2,k,ht,m^2,intmet); data(1,5) = toc;
    
    solution2 = @ (t,x,y) t*x*(x-1)*sin(pi*y);
    sol2 = zeros(m^2,k);
    for i = 1:k
        for j = 1:m
            sol2(1+(j-1)*(m):j*(m),i) = solution2(T(i),X(j+1),Y(2:end-1));
        end
    end
    Err = zeros(1,m^2);
    for i = 1:m^2
        Err(i) = norm(sol2(i,:)-U(i,:) ,Inf)/max(abs(sol2(i,:)));
    end
    data(1,4) = max(Err);
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