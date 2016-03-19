% Finite difference scheme of the diffusion equation
% -------------
% u_t = c*u_xx
% -------------
% using Crank-Nicholson and Neumann boundaries

close all;

M = 100;             % number of spacial grid points
K = 7000;           % number of time steps
h = 1/(M+1);        % space step size
k = 5e-5;           % time step size

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

N = 5000/5000;              % Scaled number of neurotransmitters
R = 152/5000;       % Scaled number og receptors, density*area = 152
RN = 0;             % Scaled number of connections

density_R = 152/5000;

epsilon = L/50;    % Extracellular fluid with
range = floor(epsilon/h);

% R + N + 2*RN = constant
% b = R + N - 2*RN; % rate of flux

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

flux = 0;
d = zeros(M+2,1);
d(M+2) = flux;

%A(M+2, M+1:M+2) = [1 -1];
A(M+2, M:M+2) = [-1/2 2 -3/2];
A(M+1,M+2) = 1;

A = scale*r/2*A;
I = eye(M+2);


% Reaction constants
% k_scale = 6.022e23/1e-3;
k1 = 1e3;
k2 = 1e1;
N_U = norm(U(1,end-range:end),1);
P_R = 1;

% figure;
figure(1)
filename = 'testnew51.gif';

for n = 1:K
    % Update flux for next time step using information from this time step:
    if U(n,end) > 0
        
        P_R = (R)/(density_R)    % Relative amount of available receptors
        if P_R < 0
            P_R = 0;
        end
        
        flux = (k1*P_R*N_U - k2*(1-P_R))*density_R
    else
        flux = 0
    end
    %d(:,n+1) = [zeros(M+1,1); flux];
    d(end-range+1:end, n+1) = flux/range;
    
    U(n+1,:) =  (I-A)\((I+A)*(U(n,:)')) - (I-A)\(k/h*(d(:,n)+d(:,n+1)));
    
    % Update number of receptors and connections:
    R = R - flux*k;
    RN = RN + flux*k;
    N_U = norm(U(n,end-range:end),1);
    
    plot([0,x,L],U(n+1,:),'.')
    hold on
    plot([x(end-range+1:end), L],0, 'r*');
    hold off
    xlabel('x')
    ylabel('Concentration of neurotransmitters')
    axis([0, L, 0, N])
    
%     drawnow
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if n == 1;
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
    
    if flux < 1e-1 && P_R < 0.8
        time = n*k
        break
    end
end



% figure;
% mesh(U)
% ylabel('time')
% xlabel('x')
% zlabel('Relative density of N')
% title('Distribution of neurotransmitters (N) over a line')
% % save('figure_1d.pdf')
