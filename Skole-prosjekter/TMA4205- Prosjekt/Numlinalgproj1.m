clear

u0=1;
un=-1;
itt = 1000;



N=20;

u = zeros(N-1,itt);

x = linspace(0,1,N+1);
h = 1/N;


alpha = -2/h^2+2/h;     
gamma = 1/h^2-2/h;  
delta = 1/h^2;
tau1  = u0*(-1/h^2+2/h);    
tau2  = -un*(1/h^2);        

A = alpha*diag(ones(N-1,1))+(delta)*diag(ones(N-2,1),1)+(gamma)*diag(ones(N-2,1),-1);

b(:,1) = -pi^2*cos(pi*x(2:N))-2*pi*sin(pi*x(2:N)); 
b(1)   = b(1)+ tau1;
b(end) = b(end)+tau2;


uex(:,1) = cos(pi*x);

uting(:,1) = A\b;  





D = alpha*diag(ones(N-1,1));


%u(:,1) = zeros(N-1,1);
u(:,1) = linspace(u0,un,N-1);
k = 2;
e(:,1) = uex-[u0;u(:,1);un];

while (k<= itt)

u(:,k) = D\(D-A)*u(:,k-1) + D\b;
e(:,k) =uex-[u0;u(:,k);un];
k = k+1;
end

figure(1)
plot(x,uex,'b')
hold on
plot(x,[u0;u(:,end);un],'r')
plot(x,[u0;uting;un],'g')

figure(2)
enorm = 1:itt;
for i = 1:itt
enorm(i) = log(norm(e(:,i),Inf));
end
plot(1:itt,enorm,'g')
hold on

