% Number of points
N = 125;



% Er det en feil i programmet) U(1) != 0 for små N.
% Fungerer greit for N >= 16



% Get grid
[p,tri,edge] = getDisk(N);
n = size(tri);

% Make zero vectors and matrices
b = zeros(N,1);
A = zeros(N,N);
U = zeros(N,1);
%Nbound = zeros(N,1);
noe = zeros(N,1);

% Pionts not on the edge
nonedge = 1:N;
nonedge(edge(:,1)) = [];
% Points with Neumann boundary
Nbound = edge(:,1);
Nbound(p(edge(:,1),2) > 0 ) = [];

% Function handlers
f = @(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2));
g = @(x,y) 4*pi*sqrt(x^2+y^2)*cos(2*pi*(x^2+y^2));
gkant = 4*pi;


display('Time it takes to make A and b')
% Assemble A and b
tic;
for i = 1:n(1)
    for k = 1:n(2)
        temp2 = coeff(tri,p,k,i);
        phi = @(x,y) temp2(1)*x+temp2(2)*y + temp2(3);
        fphi = @(x,y) f(x,y)*phi(x,y);
        b(tri(i,k)) = b(tri(i,k)) + quadrature2D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),4,fphi);
        for j = 1:n(2)
            temp1 = coeff(tri,p,j,i);
            nabsq = @(x,y) [temp2(1),temp2(2)]*[temp1(1);temp1(2)];
            A(tri(i,k),tri(i,j)) = A(tri(i,k),tri(i,j)) + quadrature2D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),1,nabsq);
        end
    end
end
toc;


display('Time it takes to put Neumann conditions in to b')
tic;
for i = edge'
    k = getTriangle(tri,i(1),i(2));     % k(1) er trekant nummer, k(2) er hjørne nummer.
    temp1 = coeff(tri,p,k(2),k(1));
    phi1 = @(x) temp1(1)*cos(x)+temp1(2)*sin(x) + temp1(3);
    gphi1 = @(x) phi1(x)*gkant;
    
    % Må transformere koordinatene til sirkelkoordinater.
    p1 = atan2(p(i(1),2),p(i(1),1));
    p2 = atan2(p(i(2),2),p(i(2),1));
    if p1 < 0
        p1 = p1 + 2*pi;
    end
    if p2 <= 0
        p2 = p2 + 2*pi;
    end
    temp1 = quadrature1D(p1,p2,4,gphi1);
    b(i(2)) = b(i(2)) + temp1*2;
end
toc;




display('Time it takes for solving A\b')
tic;
% Solving system
U([nonedge,Nbound']) = A([nonedge,Nbound'],[nonedge,Nbound'])\b([nonedge,Nbound']);
toc;

display('Time it takes for plotting')
tic;
% Element solution:
figure(1)
scatter3(p(:,1),p(:,2),U)
hold on

% Analytical solution:
figure(1)
u = @(x,y) sin(2*pi*(x.^2+y.^2));
ezmesh(u,[-1,1,-1,1])
hold off
toc;
normA = norm(A*u(p(:,1),p(:,2))-U,2)


