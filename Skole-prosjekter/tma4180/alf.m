function [x] = alf(l,x,lambda,my,g,m)
%%% Augmented lagrangian function
%%% uses the function "minimizer" to minimize a lagrangian for each
%%% iteration.

tol = 1E-6;
xold = -x +2*tol;

while norm(x-xold,2) > tol
    xold = x;
    x = minimizer(l,x,lambda,my,g,m);
    for i = 2:1:length(l)+1
        lambda(i-1) = lambda(i-1)-my*getConstrains( x,l,lambda,i );
    end
    my = 2*((my+1)^2); 
end
end


function [ c ] = getConstrains( x,l,lambda,i )
len = length(l);
y = x(len+2:end);
c = ((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2-l(i-1)^2);
end
