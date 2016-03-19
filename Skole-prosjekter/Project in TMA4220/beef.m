% Initial temperature
temperaturstart = 10;
temperaturluft = 20;
temperaturslutt = 180;
temperatursnu = 60;
temperaturferdig = temperatursnu;
snudd = 0;
itt = Inf;

% Diffusion coeffisions
alphalk = 3*10^-4;                 % Overføring mellom luft og kake
alphakk = 8.75*10^-8;                    % Overføring mellom kake og kake
alphamk = 0.7;                   % Overføring mellom metall og kake

% Number of points in space
n1 = 5;
n2 = 5;
n3 = 5;
N = n1*n2*n3;

% Number of points in time
tstart = 1;
tstep = 1;
tmax = 900;

% get grid
%[p,tri,edge] = getSphere(N);       %Nx3, 537x4, 130x3
[p,tri,edge] = getBeef(n1,n2,n3);
%p = p*1000;
n = size(tri);



% Allocating memory
A = zeros(N,N);
M = zeros(N,N);
YY = [];

display('Assembling A and M')
% Assmebling A and M	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for i = 1:n(1)
    for k = 1:n(2)
        temp2 = findPhi(p(tri(i,:),:),k);
        phi = @(x,y,z) temp2(1)*x+temp2(2)*y + temp2(3)*z + temp2(4);
        for j = 1:k
            temp1 = findPhi(p(tri(i,:),:),j);
            nphisq = @(x,y,z) [temp2(1),temp2(2),temp2(3)]*[temp1(1);temp1(2);temp1(3)];
            A(tri(i,k),tri(i,j)) = A(tri(i,k),tri(i,j)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,nphisq);
            phisq = @(x,y,z) phi(x,y,z)*(temp1(1)*x+temp1(2)*y+temp1(3)*z+temp1(4));
            M(tri(i,k),tri(i,j)) = M(tri(i,k),tri(i,j)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,phisq);
        end
    end
end
A = A+ triu(A,1)'+tril(A,-1)';
M = M+ triu(M,1)'+tril(M,-1)';
T1 = toc;
% Done assembeling A and M   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Solving M*du/dt = -alpha*u')
% Solving M*du/dt = -alpha*u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% Edge initially shared with metal
hotedge = 1:n1*n2;
% Edge initially shared with metal
coledge = unique(edge);
coledge(hotedge) = [];



% initial temperatur of the food
u0 = temperaturstart*ones(N,1);
YY(:,1) = u0;
i = 2;
%while true
for i = 2:1:floor(tmax/tstep)
    % Checking if should turn
    if u0(ceil(n1*n2*n3/2)) >= temperatursnu && ~snudd
        snudd = 1;
        ittsnu = i;
        hotedge = (N-n1*n2+1):N;
        coledge = unique((edge(1:end-n1*n2,1)));
    end
    % Change in temperature on the edges
    u0(coledge) = u0(coledge) + tstep*alphalk*(temperaturluft-u0(coledge));
    u0(hotedge) = u0(hotedge) + tstep*alphamk*(temperaturslutt-u0(hotedge));
    
    % Finding the next time step
    u0  = u0-tstep*alphakk*(M\A)*u0;
    
    % Saving the timestep
    YY(:,i) = u0;
    % Checking if done
    if snudd && min(u0) > temperatursnu
        itt = i;
        break;
    end
end

T2 = toc;
% Done Solving M*du/dt = -alpha*u%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display('plotting')
% Propper plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
% If movie
%F(tmax/tstep) = struct('cdata',[],'colormap',[]);
for i = 1:1:min(floor(tmax/tstep),itt)
    figure(1)
    scatter3( p(:,1), p(:,2), p(:,3), [], YY(:,i), 'filled' )
    title('Beef-heating'); xlabel('x');ylabel('y')
    view(-20, 40);
    caxis([min([temperaturstart, temperaturslutt,temperaturluft]),max([temperaturstart, temperaturslutt,temperaturluft])]);
    h = colorbar;
    h.Label.String = 'Temperature in {\circ} celsius';
    
    
    drawnow
    %slice(YY(:,i),p(:,1),p(:,2), p(:,3))
    % If movie
    %F(i) = getframe;
end
T3 = toc;
% If movie
%movie(F);
% done propper plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Time it takes to assemble A and M: \t\t %d \nTime it takes to solve differential equation: \t %d\nTime it takes for propper plotting: \t\t %d \n',T1,T2,T3);
fprintf('Time to turn:\t%d\nTime to done:\t%d\n',ittsnu*tstep,itt*tstep)
