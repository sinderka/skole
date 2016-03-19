clear all
close all


a = 0;
b = 8;
N = 6;

h = (b-a)/(N-1);

% mass matrix
M = 2*h/3*diag(ones(N,1)) +h/6*diag(ones(N-1,1),1) +h/6*diag(ones(N-1,1),-1);

% stiffness matrix
K = 2/h*diag(ones(N,1)) -1/h*diag(ones(N-1,1),1) -1/h*diag(ones(N-1,1),-1);

% Neumann BCs:
M(1,1) = h/3;
M(N,N) = h/3;
K(1,1) = 1/h;
K(N,N) = 1/h;

%% inhomogenous BCs:
% setting constants:
alpha = 3;  % diffusion coefficient
Pa = 1;
Pb = 1;

da = zeros(N,1);
da(1) = 1/alpha;
Qa = diag(da);

db = zeros(N,1);
db(N) = 1/alpha;
Qb = diag(db);


% initial values
U = zeros(N,1);
U(ceil(N/2)) = 1;

% building large system (incorporating boundary conditions)
Xt = [U; Pa; Pb];
Mhat = [M, zeros(N,2); zeros(2,N), [1, 0; 0, 1]];
Khat = [K, zeros(N,2); zeros(2,N+2)];
Qahat = @(X) [Qa, zeros(N,2); zeros(2,N), [X(1), 0; 0, 0]];
Qbhat = @(X) [Qb, zeros(N,2); zeros(2,N), [0, 0; 0, X(N)]];
dahat = [da; 0; 0];
dbhat = [db; 0; 1];
eahat = [zeros(N,1); 1; 0];


% settings constants
k1 = 0.5;
k2 = 0.5;
k3 = 0.5;
k4 = 0.5;
k5 = 0.5;

% function to integrate
myfun = @(t,X) -alpha*Khat*X - k3*Qahat(X)*X + k4*(1-X(N+1))*dahat ...
     +  (k4 + k5)*(1-X(N+1))*eahat - k1*Qbhat(X)*X + k2*(1-X(N+2))*dbhat;

options = odeset('Mass', Mhat);
T = 4;
[t,xa] = ode45(myfun,[0 T],Xt, options);

figure
plot (t, xa)
shading flat
legend('n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'Pa', 'Pb')
axis([0 T -0.05 1.05])
xlabel('Time')
ylabel('Concentration (Probability)')


figure
surf(a:h:b, t, xa(:,1:N))
shading flat

