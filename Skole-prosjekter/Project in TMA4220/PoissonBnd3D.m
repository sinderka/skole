% Number of points
N = 125;

% Get grid
[p ,tri,edge] = getSphere(N);
n = size(tri);

% Make zero vectors and matrices
b = zeros(N,1);
A = zeros(N,N);
U = zeros(N,1);
noe = zeros(N,1);
% Pionts not on the edge
alledge = min(min(edge)):N;
nonedge = 1:min(min(edge))-1;
nBound = find(p(alledge,1) >= 0)';

% Function handlers
f = @(x,y,z) -8*pi*cos(2*pi*(x^2+y^2+z^2))+16*pi^2*(x^2+y^2+z^2)*sin(2*pi*(x^2+y^2+z^2));

g = 4*pi;

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

%Bad way to find the b needed for the Neumann condition
display('Time it takes to put Neumann conditions in to b')
tic;
for i = edge'
    for j = 1:length(edge(1,:))
    k = getTriangle3(tri,i(1),i(2),i(3));
    temp1 = coeff3(tri,p,j,k(1));
    phi1 = @(x,y) temp1(1)*sin(x)*cos(y)+temp1(2)*sin(x)*sin(y) + temp1(3)*cos(x) + temp1(4);
    gphi1 = @(x,y) phi1(x,y)*g*sin(x);
    
    %Må transformere koordinatene til sirkelkoordinater.
    [azimuth,elevation,~] = cart2sph(p(i,1),p(i,2),p(i,3));
    azimuth = mod(azimuth+2*pi,2*pi);
    elevation = mod(elevation+pi,pi);
    
    temp1 = quadrature2D([azimuth(1),elevation(1)],[azimuth(2),elevation(2)],[azimuth(3),elevation(3)],4,gphi1);
    b(i(j)) = b(i(j)) + temp1; 
    end
end
toc;


display('Time it takes for solving A\b')
tic;
% Solving system
U([nonedge,nBound]) = A([nonedge,nBound],[nonedge,nBound])\b([nonedge,nBound]);
toc;





display('Time it takes for plotting')
tic
figure(1)
triboundary = convhull(p);
trisurf(triboundary,p(:,1),p(:,2),p(:,3),U);
figure(2)
ti = -1:0.05:1;
[qx,qy,qz] = meshgrid(ti,ti,ti);
F = TriScatteredInterp(p,U);
U_new = F(qx,qy,qz);
pat = patch(isosurface(qx,qy,qz,U_new,0.5) );
set(pat,'FaceColor','red','EdgeColor','none');
daspect([1 1 1]);
view(3);
camlight;
lighting phong

toc










