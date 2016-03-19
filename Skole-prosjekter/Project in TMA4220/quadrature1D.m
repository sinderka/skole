function [ I ] = quadrature1D(a,b,Nq,g)

l = b-a;
A = (a+b)/2;


if Nq == 1      %Sjekk 
   I=l*feval(g,A);
    
elseif Nq == 2  % Tror ok
    c1 =A+sqrt(3)/6*l;
    c2 =A-sqrt(3)/6*l;
   I = l*0.5*( feval(g,c1)+feval(g,c2));
elseif Nq == 3  % Tror ok
    c1 = A+sqrt(15)/10*l;
    c2 = A;
    c3 = A-sqrt(15)/10*l;
    
    I=l*(5/18*feval(g,c1)+4/9*feval(g,c2)+5/18*feval(g,c3));
    
    
elseif Nq == 4
    c1 = A-sqrt(525+70*sqrt(30))/70*l;
    c2 = A-sqrt(525-70*sqrt(30))/70*l;
    c3 = A+sqrt(525-70*sqrt(30))/70*l;
    c4 = A+sqrt(525+70*sqrt(30))/70*l;
    w1 = (18-sqrt(30))/72;
    w2 = (18+sqrt(30))/72;
    
    I = l*(w1*(feval(g,c1)+feval(g,c4))+w2*(feval(g,c2)+feval(g,c3)));
    
    
else
    disp('Something went wrong!')
end


end

