clear

u0=1;
un=-1;
itt = 3000;



N=40;

uJ = zeros(N-1,itt);

x = linspace(0,1,N+1);
h = 1/N;

alpha = -2/h^2+2/h;
gamma = 1/h^2-2/h;
delta = 1/h^2;
tau1  = u0*(-1/h^2+2/h);      
tau2  = -un*(1/h^2);

A = alpha*diag(ones(N-1,1))+(delta)*diag(ones(N-2,1),1)+(gamma)*diag(ones(N-2,1),-1);
E =  -(gamma)*diag(ones(N-2,1),-1);
D = alpha*diag(ones(N-1,1));
F = -(delta)*diag(ones(N-2,1),1);

b(:,1) = -pi^2*cos(pi*x(2:N))-2*pi*sin(pi*x(2:N)); 
b(1)   = b(1)+ tau1;
b(end) = b(end)+tau2;


uex(:,1) = cos(pi*x);

uting(:,1) = A\b;    

uJ(:,1) = zeros(N-1,1);
%uJ(:,1) = linspace(u0,un,N-1);
ufGS = uJ;
ubGS = uJ;

k = 2;
eJ(:,1) = uex-[u0;uJ(:,1);un];
ebGS(:,1) = uex-[u0;ubGS(:,1);un];
efGS(:,1) = uex-[u0;ufGS(:,1);un];

while (k<= itt)
uJ(:,k) = D^(-1)*(E+F)*uJ(:,k-1) + D^(-1)*b;
ufGS(:,k) = (D-E)\F*ufGS(:,k-1) + (D-E)\b;
ubGS(:,k) = (D-F)\E*ubGS(:,k-1) + (D-F)\b;


eJ(:,k) =uex-[u0;uJ(:,k);un];
efGS(:,k) = uex-[u0;ufGS(:,k);un];
ebGS(:,k) = uex-[u0;ubGS(:,k);un];

k = k+1;
end

figure(1)
plot(x,uex,'b')
hold on
plot(x,[u0;uJ(:,end);un],'r')
plot(x,[u0;ubGS(:,end);un],'k')
plot(x,[u0;ufGS(:,end);un],'c')
plot(x,[u0;uting;un],'g')


enorm = zeros(3,itt);
figure(2)
for i = 1:itt
enorm(1,i) = log(norm(eJ(:,i),Inf));
enorm(2,i) = log(norm(ebGS(:,i),Inf));
enorm(3,i) = log(norm(efGS(:,i),Inf));
end
plot(1:itt,enorm(1,:),'b')
hold on
plot(1:itt,enorm(2,:),'r')
plot(1:itt,enorm(3,:),'g')
