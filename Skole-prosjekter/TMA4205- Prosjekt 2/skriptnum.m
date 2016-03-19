%Precomputing, 
[A_x,A_y,B_x,B_y,n_x,n_y]=prep; 
alpha=1/((n_x*n_y)^2);
e=ones(n_x*n_y,1);
L_x=chol(A_x,'lower');
L_y=chol(A_y,'lower');

P=e/norm(e);

%cg in iteration
x=B_x*P;
x=forwardSubstitution(L_x,x,(n_x-1)*n_y);
x=forwardSubstitution(L_x',x,(n_x-1)*n_y);
x=B_x'*x;

y=B_y*P;
y=forwardSubstitution(L_y,y,(n_y-1)*n_x);
y=forwardSubstitution(L_y',y,(n_y-1)*n_x);
y=B_y'*y;

z=alpha*sum(P);

SP=x+y+z;


