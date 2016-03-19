function [ I ] = quadrature1D(T,p,Nq,g)

%Må skrive om slik at den tar inn noen x = n*(h)(tar inn n), samt en vector
%eller noe slikt: g

%Tar inn feks
%T = [1,2];
%p = [1,2];      
%pow = 2;
%g = [1,2^pow,3^pow];
%skal gi ut 
h = (T(1)-T(2));
A =((T(1)+T(2))/2);
%I = l*g(A);

%I = (2-1)*(g(2)-g(1))/2;
%nærmeste vi kommer er;

%l = max(p)-min(p);
%A = (max(p)+min(p))/2;


%Nq = 3;

if Nq == 1
    %%%%%I=l*g(p);           % FEIL!
    I = h*(g(p(1))+g(p(2)))/2;
    
 elseif Nq == 2
    tall = sqrt(3)/6;
    c1 =round(A+tall*h);
    c2 =round(A-tall*h);
    I = h*0.5*(g(c1)+g(c2));
 elseif Nq == 3
    tall = sqrt(15)/10;
    c1 = round(A+tall*h);
    c2 = round(A);
    c3 = round(A-tall*h);
    I=h*(5/18*g(c1)+4/9*g(c2)+5/18*g(c3));
%     
%     
% elseif Nq == 4
%     tall = sqrt(525+70*sqrt(30))/70;
%     c1 = A-tall*l;
%     c2 = A-tall*l;
%     c3 = A+tall*l;
%     c4 = A+tall*l;
%     w1 = (18-sqrt(30))/72;
%     w2 = (18+sqrt(30))/72;
%     
%     I = l*(w1*(feval(g,c1)+feval(g,c4))+w2*(feval(g,c2)+feval(g,c3)));
%     
%     
else
    disp('Something went wrong!')
end


%end

