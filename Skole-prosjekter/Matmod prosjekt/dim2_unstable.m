 % Diffusion equation on non-square scheme in 2 dimensions
 close all
 clear all

%   Lx = 0.44e-6;
%   Ly = 15e-9;

Lx = 9;
Ly = 0.3;

nx = 500;
ny = 15;

hx = Lx/(nx+1);
hy = Ly/(ny+1);

% Number of grid points including boundary points
Mx = nx + 2;
My = ny + 2;

% Time steps:
k = 1e-6;     
K = 500;

c = 150;      % Diffusion coefficient

rx = c*k/(hx^2);
ry = c*k/(hy^2);
r = c*k/(hx^2+hy^2);


N = 5000;   % Number of neurotransmitters
R = 152;    % Number of receptors
NR = 0;     % Number of bounded R-N

U = zeros(Mx*My, 1);

% N uniformly distributed on one side:
% U(end-Mx+2:end-1,1) = ones(nx,1)*N/nx;

% Put all N in one point
U(end-floor(Mx/2)) = N;

% Creating matrix system:

I_A = speye(Mx,Mx);
C = -(c*k/(hy^2))*I_A;

E = sparse(2:Mx,1:Mx-1,1,Mx,Mx);
S = (c*k/(hx^2))*(-E-E' + 4*I_A) + I_A;

% Dirichlet boundariies at x = 0 and x = Lx:
S(1,1) = 1;     S(1,2) = 0;     S(2,1) = 0;
S(Mx,Mx) = 1;   S(Mx-1,Mx) = 0; S(Mx, Mx-1) = 0;
C(1,1) = 0;     C(Mx,Mx) = 0;

Ey = sparse(2:My,1:My-1,1,My,My);

A = kron(eye(My),S) + kron(Ey+Ey',C);

% Adding Neumann boundary conditions to U(y=0) and U(y=Ly):
Neumann = 1/(2*hy)*(4*I_A - E - E');
Neumann(1,1:2) = [-1/hx, 0];   Neumann(Mx,Mx-1:Mx) = [0,-1/hx];
Neumann_2 = -2/(2*hy)*I_A;
Neumann_2(1,1) = 1/hx;     Neumann_2(Mx,Mx) = 1/hx;

A(1:Mx, 1:2*Mx) = [Neumann, Neumann_2];
A(end-Mx+1:end, end-2*Mx+1:end) = [Neumann_2, Neumann];


d = zeros(Mx*My,1);
% flux = U(1:Mx,end)+ones(Mx,1)*(R/Mx-2*NR/Mx);
% d(1:Mx,1) = flux


figure;
for n = 1:K
    d(:,n+1) = d(:,n);
    if sum(U(1:Mx,end)) > 0 && R > 0
        %flux = 1e4*norm(U(1:Mx,end))^2*(U(1:Mx,end)+ones(Mx,1)*(R/Mx-2*NR/Mx));
%         flux = 1e-1*(U(1:Mx,end)+ones(Mx,1)*(R/Mx-2*NR/Mx));
        
        flux = 1e2*norm(U(1:Mx,end))*U(1:Mx,end)*R;
        norm(flux);
        sum(U(1:Mx,end));
    	
    else
        flux = zeros(Mx,1);
    end
    d(1:Mx,n+1) = flux;
    
    
    U(:,n+1) = (A)\U(:,n) - 0.5*A\(d(:,n)+d(:,n+1));
    
    norm(U(:,n+1),1)
    
    R = R - sum(flux);
    NR = NR + sum(flux);
    
    b = sum(U(1:Mx,end)+ones(Mx,1)*(R/Mx-2*NR/Mx));
   
    U_mesh = zeros(Mx,My);
    for i = 1:Mx
        for j = 1:My
            U_mesh(i,j) = U(i+Mx*(j-1), n+1);
        end
    end
    
    mesh(U_mesh);
    xlabel('Long')
    ylabel('Short')
    title('Distribution of N on domain per time')
    pause(0.001)    
    
end
figure;
 plot(sum(d))
 title('Flux out of domain on end over time')

