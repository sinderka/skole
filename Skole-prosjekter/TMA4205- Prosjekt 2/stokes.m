function stokes
% Verification of the Stokes discretization: Kovasznay flow

% Size of the domain/mesh
Lx = 2.5;   Ly = 2.0;
nx = 25;  ny = 20;
hx = Lx/nx; hy = Ly/ny;

% analytical solution
[u,v,p,fx,fy] = kovasznay;
gx = u; gy = v;


% arrays to store the mesh size and the errors
h  = [];
eu = [];
ev = [];
ep = [];

% verification: convergence study against an analytical solution
for i=1:5,
    % setup the system
    [Ax,Ay,Bx,By] = setup_lhs(nx,ny,hx,hy);
    [Fx,Fy,G]     = setup_rhs(Lx,Ly,nx,ny,hx,hy,fx,fy,gx,gy);

    Nu = (nx-1)*ny; Nv = nx*(ny-1); 
    Np = nx*ny; N = Nu+Nv+Np;
    Iu = 1:Nu; Iv = Nu+1:Nu+Nv; Ip = Nu+Nv+1:N;

    % Full system
    A = [Ax,            sparse(Nu,Nv), Bx;       ...
         sparse(Nv,Nu), Ay,            By;       ...
         Bx.',          By.',          sparse(Np,Np)];
    F = [Fx;Fy;G];
    % Since pressures are determined only up to a constant, we can 
    % eliminate one of them fro the system to obtain a non-singular
    % matrix.
    %
    % Warning: in the project, you should follow other strategies,
    % which are more suitable to iterative solvers.
    %
    A(:,N)=0; A(N,:)=0.; A(N,N) = 1.0; F(N) = 0.0;
    % use a direct solver: in the project, other strategies should
    % be used
    %UVP = A\F;
    
    M = [Ax,            sparse(Nu,Nv), Bx;       ...
        sparse(Nv,Nu), Ay,            By;       ...
        Bx.',          By.',          speye(Np)];
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %With preconditioner
    %Axc=ichol(Ax,struct('type','ict','droptol',1e-4));Ayc=ichol(Ay,struct('type','ict','droptol',1e-4));UVP = minres(A,F,10^-10,10000, @(X)multprec(X,Axc,Ayc,Bx,By,Nu,Nv,Np )); U = UVP(Iu); V = UVP(Iv); P = UVP(Ip);
    %No preconditioner
    UVP = minres(A,F,10^-10,10000, []); U = UVP(Iu); V = UVP(Iv); P = UVP(Ip);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %U = UVP(Iu); V = UVP(Iv); P = UVP(Ip);
    
    % remember: pressure is only defined up to a constant!
    P = P - mean(P);
    
    % vizualize the solution on the coarsest grid
    if i==1,
        figure(1);
        vizualize(U,V,P,Lx,Ly,nx,ny,hx,hy,gx,gy);
    end
    
    % measure the error
    x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1);
    [Y,X] = meshgrid(y,x);
    u_exact = u(avg(X(2:end-1,:)')',avg(Y(2:end-1,:)')');
    v_exact = v(avg(X(:,2:end-1) ) ,avg(Y(:,2:end-1) ) );
    p_exact = p(avg(avg(X')'),      avg(avg(Y')'));
    % remember: pressure is only defined up to a constant!
    p_exact = p_exact(:) - mean(p_exact(:));
    h = [h,  max(hx,hy)];
    eu= [eu, norm(U-u_exact(:),inf)];
    ev= [ev, norm(V-v_exact(:),inf)];
    ep= [ep, norm(P-p_exact,   inf)];
    % increase the size of the mesh
    nx = nx*2; ny = ny*2;
    hx = hx/2; hy = hy/2;
end


% vizualize the behaviour of the error
figure(2);
loglog(h,eu,'r*-',...
       h,ev,'b*-',...
       h,ep,'k*-');
grid on;
legend('Error U','Error V','Error P', ...
       'Location', 'Best');

%----------------------------------------------------------------------
% computation of the left-hand side block matrices
%----------------------------------------------------------------------
function [Ax,Ay,Bx,By] = setup_lhs(nx,ny,hx,hy)
Ax = kron(speye(ny)*hy,K1(nx-1,hx,2))+kron(K1(ny,hy,3),speye(nx-1)*hx);
Ay = kron(speye(ny-1)*hy,K1(nx,hx,3))+kron(K1(ny-1,hy,2),speye(nx)*hx);
Bx = kron(speye(ny)*hy,K2(nx));
By = kron(K2(ny),speye(nx)*hx);
%----------------------------------------------------------------------
% LHS: auxiliary functions
%----------------------------------------------------------------------
function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 -1;ones(n-2,1)*[-1 2 -1];-1 a11 -1],-1:1,n,n)/h;
function A = K2(n)
A = spdiags(ones(n-1,1)*[-1 1],[0 1],n-1,n);
%----------------------------------------------------------------------
% computation of the right-hand side block vectors
%----------------------------------------------------------------------
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
%----------------------------------------------------------------------
% RHS, vizualization: auxiliary functions
%----------------------------------------------------------------------
function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end
%----------------------------------------------------------------------
% vizualization of the results
%----------------------------------------------------------------------
function vizualize(U,V,P,Lx,Ly,nx,ny,hx,hy,gx,gy)
x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1);
[Y,X] = meshgrid(y,x);
% evaluate the boundary term at the appropriate grid points
% {N,S,W,E} = {north, south, east, west}
uN = gx(x,Ly);      vN = gy(avg(x),Ly);
uS = gx(x,0);       vS = gy(avg(x),0);
uW = gx(0,avg(y));  vW = gy(0,y);
uE = gx(Lx,avg(y)); vE = gy(Lx,y);
% reshape the vectors
U = reshape(U,nx-1,ny);
V = reshape(V,nx,ny-1);
P = reshape(P,nx,ny);
Ue = [uS' avg([uW;U;uE]')' uN'];
Ve = [vW;avg([vS' V vN']);vE];
Len = sqrt(Ue.^2+Ve.^2+eps);
%subplot(1,3,1);
figure(1);
cla, contourf(avg(x),avg(y),P',20,'w-'), hold on;
quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-');
p = sort(P(:)); caxis(p([8 end-7])); colorbar;
axis equal, axis([0 Lx 0 Ly]); hold off;
title('P and [U,V]');
%subplot(1,3,2);
figure(3);
pcolor(x,y,Ue'); shading('interp'); colorbar;
axis equal; axis([0 Lx 0 Ly]); title('U');
figure(4);
%subplot(1,3,3);
pcolor(x,y,Ve'); shading('interp'); colorbar;
axis equal; axis([0 Lx 0 Ly]); title('V');
%----------------------------------------------------------------------
% analytical solution: Kovasznay flow
%----------------------------------------------------------------------
function [u,v,p,fx,fy]=kovasznay
lambda  = -1;
u    = @(x,y) 1-exp(lambda*(x-0.5)).*cos(2*pi*y);
v    = @(x,y) (lambda/2/pi)*exp(lambda*(x-0.5)).*sin(2*pi*y);
p    = @(x,y) 0.5*(exp(2*lambda*(x-0.5)));
d2u  = @(x,y) (-lambda^2 + 4*pi^2)*exp(lambda*(x-0.5)).*cos(2*pi*y);
d2v  = @(x,y) (lambda^3/2/pi - 2*pi*lambda)*exp(lambda*(x-0.5)).*sin(2*pi*y);
px   = @(x,y) lambda*exp(2*lambda*(x-0.5));
py   = @(x,y) zeros(size(x)).*zeros(size(y));
fx   = @(x,y) -d2u(x,y) + px(x,y);
fy   = @(x,y) -d2v(x,y) + py(x,y);

function x = multprec(X,Axc,Ayc,Bx,By,Nu,Nv,Np ) %Kan skrives smartere
U = X(1:Nu); V = X(Nu+1:Nu+Nv); P = X(end-Np+1:end);
invAxU = Axc\(Axc'\U);
invAyV = Ayc\(Ayc'\V);
BxinvAxU = Bx'*invAxU;
ByinvAyV = By'*invAyV;
u = invAxU + Axc\(Axc'\(Bx*(BxinvAxU))) + Axc\(Axc'\(Bx*(ByinvAyV))) - Axc\(Axc'\(Bx*P));
v = Ayc\(Ayc'\(By*(BxinvAxU))) + invAyV + Ayc\(Ayc'\(By*(ByinvAyV))) - Ayc\(Ayc'\(By*P));
p = -BxinvAxU-ByinvAyV+P;
x = [u;v;p];





