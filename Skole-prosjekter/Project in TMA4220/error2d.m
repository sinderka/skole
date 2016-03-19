% Number of points
N = 125;

% Get grid
[p,tri,edge] = getDisk(N);
n = size(tri);

% Make zero vectors and matrices
b1 = zeros(N,1);
b2 = zeros(N,1);
A = zeros(N,N);
U = zeros(N,1);
M = zeros(N,N);

% Points not on the edge
nonedge = 1:N;
nonedge(edge(:,1)) = [];

% Function handlers
f = @(x,y) -8*pi*cos(2*pi*(x.^2+y.^2))+16*pi^2*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2));
u = @(x,y) sin(2*pi*(x.^2+y.^2));

% Assemble A and b
tic;
for i = 1:n(1)
    for k = 1:n(2)
        temp2 = findPhi(p(tri(i,:),:),k);
        phi = @(x,y) temp2(1)*x+temp2(2)*y + temp2(3);
        fphi = @(x,y) f(x,y)*phi(x,y);
        b1(tri(i,k)) = b1(tri(i,k)) + quadrature2D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),4,fphi);
        uphi = @(x,y) u(x,y)*phi(x,y);
         b2(tri(i,k)) = b2(tri(i,k)) + quadrature2D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),4,uphi);
        for j = 1:n(2)
            temp1 = findPhi(p(tri(i,:),:),j);
            nabsq = @(x,y) [temp2(1),temp2(2)]*[temp1(1);temp1(2)];
            A(tri(i,k),tri(i,j)) = A(tri(i,k),tri(i,j)) + quadrature2D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),1,nabsq);
            phisq = @(x,y) phi(x,y)*(temp1(1)*x+temp1(2)*y+temp1(3));
            M(tri(i,k),tri(i,j)) = M(tri(i,k),tri(i,j)) + quadrature2D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),4,phisq);
        end
    end
end
T1= toc;

uu = u(p(:,1),p(:,2));


figure(1)
scatter(p(:,1),p(:,2),[],A*uu-b1,'filled' );
title('Error plot for A'); xlabel('x');ylabel('y');zlabel('z');
h = colorbar;
h.Label.String = 'Au-b_1';
normA = norm(A*uu-b1,2);


figure(2)
scatter(p(:,1),p(:,2),[],M*uu-b2,'filled' );
title('Error plot for M'); xlabel('x');ylabel('y');zlabel('z');
h = colorbar;
h.Label.String = 'Mu-b_2';
normM = norm(M*uu-b2,2);

fprintf('Time it takes to assemble A, M, b_1 and b_2: \t %d \nE_A: \t\t\t\t\t\t %d\nE_M: \t\t\t\t\t\t %d \n ',T1,normA,normM);



