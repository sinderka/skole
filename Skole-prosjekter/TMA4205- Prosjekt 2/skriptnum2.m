[A_x,A_y,B_x,B_y,n_x,n_y]=prep;
hold off
vec=[20,25,30,35,40];
vec2=['r','b','y','g','k'];
vec3=zeros(5);
x=[10^-5,10^-3,10^-1,10^1,10^3];
for i=1:5
    for j=1:5
Lx = 2.5;   Ly = 2.0;
nx = vec(j);  ny = vec(i);
hx = Lx/nx; hy = Ly/ny;
[Ax,Ay,Bx,By] = setup_lhs(nx,ny,hx,hy);


alpha=1/((n_x*n_y)^2);
e=ones(n_x*n_y,1);
L_x=chol(A_x,'lower');
L_y=chol(A_y,'lower');
temp=L_x\B_x;
temp2=L_y\B_y;
temp=(L_x')\temp;
temp2=(L_y')\temp2;
S=B_x'*temp+B_y'*temp2;
kond=max(abs(eigs(S+alpha*ones(n_x*n_y))))/min(abs(eigs(S+alpha*ones(n_x*n_y))));
vec3(i,j)=kond;
    end
    figure
    loglog(x,vec3(i,:),vec2(i));
end
max(abs(eigs(S+alpha*ones(n_x*n_y))))
min(abs(eigs(S+alpha*ones(n_x*n_y))))