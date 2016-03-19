function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 -1;ones(n-2,1)*[-1 2 -1];-1 a11 -1],-1:1,n,n)/h;
end