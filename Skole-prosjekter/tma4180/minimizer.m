function [ x ] = minimizer(l,x,lambda,my,g,m)
%%% minimizes the function made in funck(...), with hessian and gradient as
%%% in hessian(...) and nab_vec(...), respectivly.
%%% Uses "Newtons method", and "steepest descent method", where newton fails.
%%% Also uses backtracking.
tol = 1E-6;
rho = 0.5;
c = 0.001;

len = length(l);
xold = -x +2*tol;
%nitt = 3000; i = 0;
while norm(xold-x,2)>tol
    
    
    

    Nabla = nab_vec(m,x,lambda,l,my,g,len);
    Hessian = hessian(m,x,lambda,l,my,g,len);
    gradient =Nabla([2:len,len+3:end-1]);
    hess = Hessian([2:len,len+3:end-1],[2:len,len+3:end-1]);
    
    

    p = zeros(length(x),1);
    p([2:len,len+3:end-1]) = -hess\gradient;
    
    if p([2:len,len+3:end-1])'*gradient >= 0
        p([2:len,len+3:end-1]) = -gradient;
    end
    
    alpha = 1;
    while funck(x + alpha*p,l,lambda,g,m,my) > funck(x,l,lambda,g,m,my) + c*alpha*gradient'*p([2:len,len+3:end-1])
       alpha = rho*alpha;
    end
    

     xold = x;
    x([2:len,len+3:end-1]) = x([2:len,len+3:end-1]) +alpha*p([2:len,len+3:end-1]);
    %i = i+1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient of the augmented lagrangian function
function nab = nab_vec(m,x,lambda,l,my,g,len)
nab = zeros(2*len+2,1);
for i = 2:1:len
    nab(i) =  derx(m,x,lambda,l,my,i,g,len);
    nab(i+len+1)= dery(m,x,lambda,l,my,i,g,len);
end
end

function s = derx(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = -lambda(i-1)*2*(x(i)-x(i-1)) ...
    -lambda(i)*2*(x(i+1)-x(i))*(-1) ...
    +my*((x(i)-x(i-1))^2+(y(i)-y(i-1))^2-l(i-1)^2)*(2*(x(i)-x(i-1))) ...
    +my*((x(i+1)-x(i))^2+(y(i+1)-y(i))^2-l(i)^2)*(2*(x(i+1)-x(i)))*(-1);
end

function s = dery(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s= m(i-1)*g*0.5+m(i)*g*0.5 ...
    -lambda(i-1)*2*(y(i)-y(i-1)) ...
    -lambda(i)*2*(y(i+1)-y(i))*(-1) ...
    +my*((x(i)-x(i-1))^2+(y(i)-y(i-1))^2-l(i-1)^2)*(2*(y(i)-y(i-1))) ...
    +my*((x(i+1)-x(i))^2+(y(i+1)-y(i))^2-l(i)^2)*(2*(y(i+1)-y(i)))*(-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hessian of the augmented lagrangian function
function Hessian = hessian(m,x,lambda,l,my,g,len)
Hessian = zeros(2*len+2);

for i = 2:1:len
    
    % Diagonal
    Hessian(i,i) = hessdiagx(m,x,lambda,l,my,i,g,len);
    Hessian(i+len+1,i+len+1) = hessdiagy(m,x,lambda,l,my,i,g,len);
    % Diagonal
    % Nedre diagonal
    if i ~= len
        Hessian(i+1,i) = hessianndiagx(m,x,lambda,l,my,i,g,len);
        Hessian(i+len+2,i+len+1) = hessianndiagy(m,x,lambda,l,my,i,g,len);
    end
    
    if i  ~= len
        Hessian(i,i+1) = hessianndiagx(m,x,lambda,l,my,i,g,len);
        Hessian(i+len+1,i+len+2) = hessianndiagy(m,x,lambda,l,my,i,g,len);
    end
    
    % Kryssledd
    Hessian(i,i+len+1) = hessianxy(m,x,lambda,l,my,i,g,len);
    Hessian(i+len+1,i) = hessianxy(m,x,lambda,l,my,i,g,len);
    if i ~= 2
        Hessian(i,i+len) = hessianx1y(m,x,lambda,l,my,i,g,len);
        Hessian(i+len,i) = hessianx1y(m,x,lambda,l,my,i,g,len);
    end
    
    if i ~= len
        Hessian(i,i+len+2) = hessianxy1(m,x,lambda,l,my,i,g,len);
        Hessian(i+len+2,i) = hessianxy1(m,x,lambda,l,my,i,g,len);
    end
    
    % Kryssledd
    
    
    
end
end

%%%


%%% Hessian diagonalen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = hessdiagx(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = -2*lambda(i-1) - 2*lambda(i) ...
    +my*((x(i)-x(i-1))^2+(y(i)-y(i-1))^2-l(i-1)^2)*2 ...
    +my*4*(x(i)-x(i-1))^2 ...
    +my*((x(i+1)-x(i))^2+(y(i+1)-y(i))^2-l(i)^2)*2 ...
    +my*4*(x(i+1)-x(i))^2;
end

function s = hessdiagy(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = -2*lambda(i-1) - 2*lambda(i) ...
    +my*((x(i)-x(i-1))^2+(y(i)-y(i-1))^2-l(i-1)^2)*2 ...
    +my*4*(y(i)-y(i-1))^2 ...
    +my*((x(i+1)-x(i))^2+(y(i+1)-y(i))^2-l(i)^2)*2 ...
    +my*4*(y(i+1)-y(i))^2;
end
%%% Hessian diagonalen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Hessian Nedre diagonal
function s = hessianndiagx(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = +2*lambda(i-1) ...
    +my*((x(i)-x(i-1))^2+(y(i)-y(i-1))^2 - l(i-1)^2)*2*(-1) ...
    +my*4*(x(i)-x(i-1))^2*(-1);
end

function s = hessianndiagy(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = +2*lambda(i-1) ...
    +my*((x(i)-x(i-1))^2+(y(i)-y(i-1))^2 - l(i-1)^2)*2*(-1) ...
    +my*4*(y(i)-y(i-1))^2*(-1);
end


%%% Hessian Kryssledd
function s = hessianxy(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = my*2*(y(i)-y(i-1))*2*(x(i)-x(i-1))+my*2*(y(i)-y(i+1))*2*(x(i)-x(i+1));
end
function s = hessianx1y(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = my*2*(y(i)-y(i-1))*(-1)*2*(x(i)-x(i-1));
end

function s = hessianxy1(m,x,lambda,l,my,i,g,len)
y = x(len+2:end);
s = my*2*(y(i+1)-y(i))*(-1)*2*(x(i+1)-x(i));
end
%%% Hessian Kryssledd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% augmented lagrangian function
function tall = funck(x,l,lambda,g,m,my)
tall = 0;
for i = 2:1: length(lambda)+1
    tall = tall +getFunc( l,x,g,m,i ) - lambda(i-1)*getConstrains( x,l,i ) +my/2*(getConstrains( x,l,i ))^2;
end


end

function [ c ] = getConstrains( x,l,i )
len = length(l);
y = x(len+2:end);
c = ((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2-l(i-1)^2);
end


function [ f ] = getFunc( l,x,g,m,i )
len = length(l);
y = x(len+2:end);
f = m(i-1)*g*(y(i)+y(i-1))/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
