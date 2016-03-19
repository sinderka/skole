% Gjør dimensjonsløst

tic;
% Step size % try to keep ht/hx^2 < 0.5
nx=10;                      % 10
ny=10;                      % 10
%nz = 10;
nt=400000;                  % 400000


my = 10^-6;

% Simulation length
time=1;

% Size of domain
d=2*0.22*my;
%h= 1;                   % Skal være 10^-9;

% Step length
hx = d/nx;
hy = d/ny;
%hz = h/nz;
ht = time/nt;


% Constants 
kappa=0.3*(my)^2;           % Bør sikkert være mindre enn 1!
k1 = 10000;                     % har ingen anelse hva disse skal være
k2 = 10;                     % har ingen anelse hva disse skal være
numR = 1000*d^2*my^(-2);    % Number of receptors
first = 1;                  % 
numN = 5000;                % Number of Nurotransmitters
TIME = 0;

dx = 1/((hx)^2);
dy = 1/((hy)^2);
%dz = 1/((hz)^2);

disk = kappa*ht*max(dx,dy); % Denne burde ikke være for stor (<1/2)

% Making A
T=((-2*dx-2*dy)*diag(ones(nx,1)) +dx*diag(ones(nx-1,1),1) +dx*diag(ones(nx-1,1),-1));
A=blktridiag(T,+dy*eye(nx),+dy*eye(nx),ny);
%A = blktridiag(A,+dz*eye(ny),+dz*eye(ny),nz);

% Allocating memory
U=zeros(nx*ny,nt);                      % For [N]
U(ceil((nx*ny)/2+nx/2),1)=numN;        % Putting all in ~the middle
%U(randi([1,nx*ny],1,1),1)=numN;        % Putting all in a random place
%U(ceil((nx*ny)/2),1)=numN;              % Putting all on the edge
R = zeros(nx*ny,nt);                    %For [R]  
R(:,1) = numR/(nx*ny)*ones(nx*ny,1);    % Start with lots of free receptors
T1 = toc;
tic;

for i=1:nt
    % Solving the differential equation for [N]
    temp2 = U(:,i)+(ht/2)*(kappa*A*U(:,i)-k1*R(:,i).*U(:,i)+k2*(numR/(nx*ny)-R(:,i)));
    pos = find(temp2>0);
    neg = find(temp2<0);
    U(pos,i+1) = temp2(pos);
    U(neg,i+1) = 0;
    % Solving the differential equation for [R]
    temp1 =R(:,i) + ht*(-k1*R(:,i).*U(:,i)+k2*(numR/(nx*ny)-R(:,i)));
    pos = find(temp1>0);
    neg = find(temp1<0);
    R(pos,i+1) = temp1(pos);
    R(neg,i+1) = 0;
    if  (sum(R(:,i+1)) > sum(R(:,i))) && first
        TIME = i*ht;
        first = 0;      
        break;
    end
    
end
T2= toc;
tic;
% For plotting over the correct domain
x = linspace(0,d,nx);
y = linspace(0,d,ny);

% Allocating memmory for plots
u = zeros(nx,ny);
RR= zeros(nx,ny);

% No need to plot all timesteps, increase nF!
nF = 1000;
for j = 1:nF:min(nt,round(TIME/ht))
    
    for i = 1:ny
         u(:,i) = U((i-1)*nx+1:i*nx,j);
        RR(:,i) = R((i-1)*nx+1:i*nx,j);
    end
    figure(1);
    surf(y,x,u); title('[N]'); %axis([0,d,0,d,0,max(max(u))]);
	drawnow
    figure(2)
    surf(x,y,RR); title('[R]'); %axis([0,d,0,d,0,numR/(nx*ny)]);
	drawnow


end
T3 = toc;

fprintf('Time init: %d\nTime loop: %d\nTime plot: %d\nTime sign: %d\n',T1,T2,T3,TIME)



