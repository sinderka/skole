function [Zn] = intmeth(F,H,hn,k,ht,m,intmet)
%%%%%Solves the equation z'-H*z=hn*F Numerically
%%% ht is stepsize, m is the size of the square matrix H
%%% len_T is the number of time-steps, zn is a vector of length len_T
%%% intmet says something about which integration-method that should be
%%% used. Boundary equal to zero.

if intmet == 1
    Zn = zeros(m,k);
    for j = 2:1:k
        Zn(:,j) = Zn(:,j-1) + H*ht*Zn(:,j-1);
        Zn(1,j) = Zn(1,j) + ht*hn*F(j-1);
    end
elseif intmet == 2
    Zn = zeros(m,k);
    Zn(1,2) = ht/2*(F(1)*hn+F(2)*hn);
    for j = 3:1:k
        Zn(:,j) = Zn(:,j-1)+ht/2*H*Zn(:,j-1);
        Zn(1,j) = Zn(1,j) + ht/2*(F(j-1)*hn+F(j)*hn);
    end
else
    disp('error in "intmeth"');
end
end