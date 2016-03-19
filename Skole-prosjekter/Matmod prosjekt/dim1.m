% This script models a neurotransmission in the brain in 1 dimension 
% using a finite difference scheme of the diffusion equation in 1d
% and and Neumann boundaries
% -------------------------------------------------------------------
% The resulting plots:
% - Distribution of neurotransmitters between axon terminal and dendritic
%   spine over time.
% - Scaled amount of flux through the end (dendritic spine) per time.
% - A mesh of the distribution with space and time axes.
% -------------------------------------------------------------------

close all;

M = 100;            % number of spacial grid points
K = 12000;          % number of time steps
h = 1/(M+1);        % space step size
k = 1e-4;           % time step size

r = k/(h^2);
L = 1;

% Grid points in space, not including boundary points
x = h:h:(1-h);

% Total time
time = K*k;        


% Scaling gives u_t = scale*u_xx
% x* = Lx, t* = Tt, u* = Nu
% Assume L = 15e-9, T = 10e-3, N = 5000
% scale = diff_const*T/L^2
scale = 4/3;

% Reaction constants
k1 = 1e1;
k2 = 1e0;

% Initial values:
N = 5000/5000;      % Scaled number of neurotransmitters
R = 152/5000;       % Scaled number og receptors, density*area = 152
RN = 0;             % Scaled number of connections

density_R = 152/5000;

% Initial value, start with all neurotransmitters in x=0:
U = zeros(1,M+2);
U(1) = N;

% Extracellular fluid width
epsilon = L/100;    
range = floor(epsilon/h);
if range < 1
    range = 1;
end


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
flux = 0;
d = zeros(M+2,1);
d(M+2) = flux;

A(M+2, M:M+2) = [-1/2 2 -3/2];
A(M+1,M+2) = 1;


A = scale*r/2*A;
I = eye(M+2);

N_U = norm(U(1,end-range:end),1);
P_R = (R)/(density_R);

figure;

for n = 1:K
    
    % Update flux for next time step using information from this time step:
    if U(n,end) > 0
        P_R = (R)/(density_R);    % Relative amount of available receptors
        if P_R < 0
            P_R = 0;
        end
        flux = (k1*P_R*N_U - k2*(1-P_R))*density_R;
        flux_t = flux*k;    % Flux per time step
    else
        flux = 0;
    end
    
    % Update Neumann boundary vector:
    d(end-range+1:end, n+1) = flux/range;
    
    % Update distribution of neurotransmitters:
    U(n+1,:) =  (I-A)\((I+A)*(U(n,:)')) - (I-A)\(k/h*(d(:,n)+d(:,n+1)));
    
    % Update number of receptors and connections:
    R = R - flux*k;
    RN = RN + flux*k;
    N_U = norm(U(n,end-range:end),1);
    
    
    if mod(n,20) == 0
        
        subplot(2,1,1)
        plot([0,x,L],U(n+1,:),'.')
        title('Distribution of N on a line over time')
        hold on
        plot([x(end-range+1:end), L],0, 'r*');
        hold off
        xlabel('z')
        ylabel('Concentration of N')
        axis([0, L, 0, N])

        subplot(2,1,2)
        plot(sum(d)*10)
        title('Scaled flux through end per time')
        axis([0, K, 0, 1])
        xlabel('Iteration')
        ylabel('Scaled flux')
        pause(0.00001)
    end

    % Register time of equilibrium:
    if flux*k < 1e-8 && P_R < 0.8 && time == k*K
        time = n*k;
    end
end
time_t = 0:k:k*K;
z = 0:h:L;

figure;
mesh(z, time_t ,U)
ylabel('time')
xlabel('z')
zlabel('Relative density of neurotransmitters')

disp('Time before signal is sent in ms: ')
disp(time)
disp('\nAmount of receptors occupied: ')
disp(1-P_R)
