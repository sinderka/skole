clear all
close all

%% Doing 2d case:
xdir = [0, 1];
ydir = [0, 1];

N1 = 10;    % x-dir
N2 = 10;    % y-dir
N = N1*N2;
Nn = (N1-2)*(N2-1);

h = (xdir(2) - xdir(1))/(N1 - 1);
l = (ydir(2) - ydir(1))/(N2 - 1);

% making points:
z = zeros(Nn, 2);

% building quads of points. Each quad is one element.
quad = zeros(1,4);
count = 0;
for j = 1:N2-1
    for i = 1:N1-2
        count = count +1;
        if count + N1 - 1 <= Nn
            quad(count,:) = [count, count + 1, count + N1-1, count + N1-2];
            %quad(count,:) = [count + N1-2, count + N1-1, count, count + 1];
        end
        z(count,:) = [(i-1)*h, (j-1)*l];
    end
end
% removing uneccesary quads:
quad = quad(mod(quad(:,1),N1-2) ~= 0,:);

% Looping over elements to find mass and stiffness matrix
M = zeros(Nn);
K = zeros(Nn);
for i = 1:length(quad)
    nodes = quad(i,:);

    Me = diag(h*l/9*ones(4,1)) ...
        + diag(h*l/18*ones(3,1),-1) + diag(h*l/18*ones(3,1),+1) ...
        + diag(h*l/36*ones(2,1), -2) + diag(h*l/36*ones(2,1), 2) ...
        + diag(h*l/18, -3) + diag(h*l/18,3);
    
    
    M(nodes, nodes) = M(nodes, nodes) + Me;
    
    Ke = 0.5*(l/(3*h) + h/(3*l))*diag(ones(4,1));  % 0.5 so wont get factor of 2 when Ke = Ke + Ke'
    Ke(1,2) = h/(6*l) - l/(3*h);
    Ke(2,3) = Ke(1,2);
    Ke(3,4) = l/(6*h) - h/(3*l);
    Ke(1,4) = Ke(3,4);
    Ke(1,3) = -h/(6*l) - l/(6*h);
    Ke(2,4) = Ke(1,3);

    Qe = Ke + Ke';
    
    K(nodes, nodes) = K(nodes, nodes) + Qe;
end

% Boundary conditions for a-side:
for i = 1:N1
    M(i,i) = 2*h*l/9;
    K(i,i) = 2*(l/(3*h) + h/(3*l));
end
% dirichlet conditions are built into system


% setting up initial values
U = [ones(N1-2,1); zeros(Nn-N1+2,1)];
PR = ones(N1,1);
Xt = [U ; PR];


% Building large system (building in all boundary conditions)
Mbar = diag(h*ones(N1,1));
Mbar(1,1) = h/2;
Mbar(N1,N1) = h/2;

alpha = 1;
% put in gammas in c
c = h/alpha*(1 + 6 + 1)*ones(Nn,1);
c(1) = h/alpha*4;
c(Nn) = h/alpha*4;      % not sure if we should be shorter?

% dimension of all bars = N1
cbar = h*ones(N1,1);
cbar(1) = h/2;
cbar(N1-1) = h/2;

% building Qbar and Q
Qbar = @(U) h/8*diag( [U(2), 6*U(2) + U(3), (U(2:N1-3) + 6*U(3:N1-2) + U(4:N1-1))', 6*U(N1-1) + U(N1-2), U(N1-1)]);
Q = @(PR) h/(24*alpha)*diag([ PR(1:N1-2) + 14*PR(2:N1-1) + PR(3:N1); zeros(Nn-N1+2,1)]) ...
    + h/(12*alpha)*diag([ PR(1:N1-3) + PR(2:N1-2); zeros(Nn-N1+2,1)],-1) ...
    + h/(12*alpha)*diag([ PR(2:N1-2) + PR(3:N1-1); zeros(Nn-N1+2,1)],1);

% d-vector:
d = @(PR) h/8*[ PR(1:N1-2) + 6*PR(2:N1-1) + PR(3:N1); zeros(Nn-N1+2,1)];

% hat members are used in final computation
chat = [c; cbar];
Mhat = [M, zeros(Nn, N1); zeros(N1, Nn), Mbar];
Khat = [K, zeros(Nn, N1); zeros(N1, Nn), zeros(N1, N1)];
Dhat = [zeros(Nn, Nn), zeros(Nn,N1); zeros(N1, Nn), Mbar];
Qhat = @(X) [Q(X(Nn+1:end)), zeros(Nn, N1); zeros(N1,Nn), Qbar(X(1:Nn))];
dhat = @(X) [d(X(Nn+1:end)); zeros(N1,1)];


% setting constants:
k1 = 0.5;
k2 = 0.5;

% solving system:
myfun = @(t,X) -alpha*Khat*X - k1*Qhat(X)*X + k2*chat - k2*dhat(X) - k2*Dhat*X;

options = odeset('Mass', Mhat);

[t,xa] = ode45(myfun,[0 0.3],Xt, options);


%% Plotting results
results = [z(:,1), z(:,2), xa(end,1:Nn)'];

figure
subplot(1,2,1)
tetramesh(quad, results)
title('solution at time 0.3')
axis square
xlabel('x'), ylabel('y'),
subplot(1,2,2)
tetramesh(quad, [z, Xt(1:Nn)])
title('initial condition')
axis square
xlabel('x'), ylabel('y'),

