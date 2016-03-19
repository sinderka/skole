% Finite difference scheme of the diffusion equation
% -------------
% u_t = c*u_xx
% -------------
% using Crank-Nicholson and Neumann boundaries

close all;

M = 80;             % number of spacial grid points
K = 5000;           % number of time steps
h = 1/(M+1);        % space step size
k = 1e-6;           % time step size

r = k/(h^2);
L = 1;



% Scaling gives u_t = scale*u_xx
% x* = Lx, t* = Tt, u* = Nu
% Assume L = 15e-9, T = 10e-4, N = 5000
% scale = diff_const*N*T/L^2
T = 10e-4;
scale = 2000/3;

% Grid points in space, not including boundary points
x = h:h:(1-h);


% -------------
% Flux
% -------------

N = 1;              % Scaled number of neurotransmitters
R = 152/5000;       % Scaled number og receptors, density*area = 152
RN = 0;             % Scaled number of connections

% R + N + 2*RN = constant
b = R + N - 2*RN; % rate of flux



% Initial value, start with all neurotransmitters in x=0:
U = zeros(1,M+2);
U(1) = N;


% --------------------
% Making matrix system:
% --------------------

ee = ones(M, 1);
e1 = ones(M-1, 1);

% Matrix for internal nodes
A_internal = sparse(diag(e1,-1) - 2*diag(ee) + diag(e1,1));

% Matrix for whole system including boundary points
A = spalloc(M+2,M+2, 3*(M+2));
A(2:M+1,2:M+1) = A_internal;

% Homogenous Neumann boundary conditions in x=0:
A(1, 1:3) = [-3/2*h, 2*h, -h/2];
A(2,1) = 1;

% Neumann boundary conditions in x=L describing rate of reaction:
% Flux = amount/t* = (1/T)*amount/t
flux = b/T*(U(1,M+2) + R - 2*RN);

d = zeros(M+2,1);
d(M+2) = flux;

A(M+2, M+1:M+2) = [2 -2];
A(M+1,M+2) = 1;

A = scale*r/2*A;
I = eye(M+2);




figure;

for n = 1:K
    % Update flux for next time step using information from this time step:
    flux = b/T*(U(n,M+2) + R - 2*RN);
    d(:,n+1) = [zeros(M+1,1); flux];
    
    % Update number of receptors and connections:
    R = R - flux*k;
    RN = RN + flux*k;
    
    U(n+1,:) =  (I-A)\((I+A)*(U(n,:)')) - (I-A)\(k/h*(d(:,n)+d(:,n+1)));
    
    
    
    plot([0,x,L],U(n+1,:),'*')
    xlabel('x')
    ylabel('Concentration of neurotransmitters')
    axis([0, L, 0, N])
    pause(0.001)
    
end

figure;
mesh(U)
ylabel('time')
xlabel('x')

