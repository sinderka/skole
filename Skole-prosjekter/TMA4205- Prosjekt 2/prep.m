function [Ax,Ay,Bx,By,nx,ny]=prep
Lx = 2.5;   Ly = 2.0;
nx = 25;  ny = 20;
hx = Lx/nx; hy = Ly/ny;
[Ax,Ay,Bx,By] = setup_lhs(nx,ny,hx,hy);
end