function [Fx,Fy,G] = setup_rhs(Lx,Ly,nx,ny,hx,hy,fx,fy,gx,gy)
x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1);
[Y,X] = meshgrid(y,x);
% evaluate the forcing term at the appropriate grid points
Fx  = [fx(avg(X(2:end-1,:)')',avg(Y(2:end-1,:)')')]*hx*hy;
Fy  = [fy(avg(X(:,2:end-1) ) ,avg(Y(:,2:end-1) ) )]*hx*hy;
% evaluate the boundary term at the appropriate grid points
% {N,S,W,E} = {north, south, east, west}
uN = gx(x,Ly);      vN = gy(avg(x),Ly);
uS = gx(x,0);       vS = gy(avg(x),0);
uW = gx(0,avg(y));  vW = gy(0,y);
uE = gx(Lx,avg(y)); vE = gy(Lx,y);
% pack the boundary contributions
Ubc = ([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']*hx/hy+...
       [uW;zeros(nx-3,ny);uE]*hy/hx);
Vbc = ([vS' zeros(nx,ny-3) vN']*hx/hy+...
       [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]*hy/hx);
Pbc = ([-uW;zeros(nx-2,ny);uE]*hy+...
       [-vS' zeros(nx,ny-2)  vN']*hx);
% combine the boundary and right-hand side contributions
Fx = Fx(:) + Ubc(:);
Fy = Fy(:) + Vbc(:);
G  = Pbc(:);
end