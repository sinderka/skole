% Number of points
N = 125;

% Get grid
[p,tri,edge] = getSphere(N);
n = size(tri);

% Make zero vectors and matrices
b1 = zeros(N,1);
A = zeros(N,N);
U = zeros(N,1);
M = zeros(N,N);
b2 = zeros(N,1);

% Pionts not on the edge
nonedge = 1:N;
nonedge([edge(:,1);edge(:,2);edge(:,3)]) = [];

% Function handlers
f = @(x,y,z) -8*pi*cos(2*pi*(x.^2+y.^2+z.^2))+16*pi^2*(x.^2+y.^2+z.^2).*sin(2*pi*(x.^2+y.^2+z.^2));
u = @(x,y,z) sin(2*pi*(x.^2+y.^2+z.^2));

%Assemble A and b
tic;
for i = 1:n(1)
    for k = 1:n(2)
        temp2 = findPhi(p(tri(i,:),:),k);
        phi = @(x,y,z) temp2(1)*x+temp2(2)*y + temp2(3)*z + temp2(4);
        fphi = @(x,y,z) f(x,y,z)*phi(x,y,z);
        uphi = @(x,y,z) u(x,y,z)*phi(x,y,z);
        b1(tri(i,k)) = b1(tri(i,k)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,fphi);
        b2(tri(i,k)) = b2(tri(i,k)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,uphi);
        for j = 1:k
            temp1 = findPhi(p(tri(i,:),:),j);
            nphisq = @(x,y,z) [temp2(1),temp2(2),temp2(3)]*[temp1(1);temp1(2);temp1(3)];
            A(tri(i,k),tri(i,j)) = A(tri(i,k),tri(i,j)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),1,nphisq);
            phisq = @(x,y,z) phi(x,y,z) *(temp1(1)*x+temp1(2)*y+temp1(3)*z+temp1(4));
            M(tri(i,k),tri(i,j)) = M(tri(i,k),tri(i,j)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,phisq);
        end
    end
end
A = A+ triu(A,1)'+tril(A,-1)';
M = M+ triu(M,1)'+tril(M,-1)';
T1 = toc;



uu = u(p(:,1),p(:,2),p(:,3));

figure(3)
scatter3(p(:,1),p(:,2),p(:,3),[],A*uu-b1,'filled' );
title('Error plot for A'); xlabel('x');ylabel('y');zlabel('z');
h = colorbar;
h.Label.String = 'Au-b_1';
normA = norm(A*uu-b1,2);


figure(4)
scatter3(p(:,1),p(:,2),p(:,3),[],M*uu-b2,'filled' );
title('Error plot for M'); xlabel('x');ylabel('y');zlabel('z');
h = colorbar;
h.Label.String = 'Mu-b_2';
normM = norm(M*uu-b2,2);

%T2 =toc

fprintf('Time it takes to assemble A, M, b_1 and b_2: \t %d \nE_A: \t\t\t\t\t\t %d\nE_M: \t\t\t\t\t\t %d \n ',T1,normA,normM);

