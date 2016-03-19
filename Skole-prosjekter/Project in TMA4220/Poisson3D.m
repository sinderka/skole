
% Number of points
N = 125;

% Get grid
[p,tri,edge] = getSphere(N);
n = size(tri);

% Make zero vectors and matrices
b = zeros(N,1);
A = zeros(N,N);
U = zeros(N,1);

% Pionts not on the edge
nonedge = 1:N;
nonedge(edge(:,1)) = [];

% Function handlers
f = @(x,y,z) -12*pi*cos(2*pi*(x^2+y^2+z^2))+16*pi^2*(x^2+y^2+z^2)*sin(2*pi*(x^2+y^2+z^2));  % hmm, skal det være pi^2 eller ikke?


display('Time it takes to make A and b')
% Assemble A and b
tic;
for i = 1:n(1)
    for k = 1:n(2)
        temp2 = findPhi(p(tri(i,:),:),k);
        phi = @(x,y,z) temp2(1)*x+temp2(2)*y + temp2(3)*z + temp2(4);
        fphi = @(x,y,z) f(x,y,z)*phi(x,y,z);
        b(tri(i,k)) = b(tri(i,k)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,fphi);
        for j = 1:k
            temp1 = findPhi(p(tri(i,:),:),j);
            nphisq = @(x,y,z) [temp2(1),temp2(2),temp2(3)]*[temp1(1);temp1(2);temp1(3)];
            A(tri(i,k),tri(i,j)) = A(tri(i,k),tri(i,j)) + quadrature3D(p(tri(i,1),:),p(tri(i,2),:),p(tri(i,3),:),p(tri(i,4),:),5,nphisq);
        end
    end
end
A = A+ triu(A,1)'+tril(A,-1)';
toc;




display('Time it takes for solving A\b')
tic;
% Solving system
U(nonedge) = A(nonedge,nonedge)\b(nonedge);
toc;




display('Time it takes for plotting')
tic
figure(1)
triboundary = convhull(p);
trisurf(triboundary,p(:,1),p(:,2),p(:,3),U);
figure(3)
u = @(x,y,z) sin(2*pi*(x.^2+y.^2+z.^2));
U = u(p(:,1),p(:,2),p(:,3));
ti = -1:0.05:1;
[qx,qy,qz] = meshgrid(ti,ti,ti);
F = TriScatteredInterp(p,U);
U_new = F(qx,qy,qz);
pat = patch(isosurface(qx,qy,qz,U_new,0.5) );
set(pat,'FaceColor','red','EdgeColor','none');
daspect([1 1 1]);
view(3);
%axis([-1,1,0,1,-1,1])
camlight;
lighting phong
toc









