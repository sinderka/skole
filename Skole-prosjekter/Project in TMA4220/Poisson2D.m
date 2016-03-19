% Number of points
N = 125;

[p tri edge] = getDisk(N);
n = size(tri);

% Make zero vectors and matrices
b = zeros(N,1);
A = zeros(N,N);
U = zeros(N,1);

% Points not on the edge
nonedge = 1:N;
nonedge(edge(:,1)) = [];

% Function handlers
f = @(x,y) -8*pi*cos(2*pi*(x.^2+y.^2))+16*pi^2*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2));
u = @(x,y) sin(2*pi*(x.^2+y.^2));

% Assemble A and b
display('Time it takes to make A and b')
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

% Solving system
display('Time it takes for solving A\b')
tic;

U(nonedge) = A(nonedge,nonedge)\b(nonedge);
toc;



display('Time it takes for plotting')
tic;
% Element solution:
figure(1)
scatter3(p(:,1),p(:,2),U); xlabel('x');ylabel('y');zlabel('z')
hold on

norm(u(p(:,1),p(:,2))-U,2)

% Analytical solution:
figure(1)
ezmesh(u,[-1,1,-1,1])
hold off
toc;



