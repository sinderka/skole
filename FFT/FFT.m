% Number of points:
N = 8000; % Number of points in func (for recording: 22050)
K = 400;% Number of coeff in estimate (for recording: 400)

if K> N/2
    K = floor(N/2);
end

% Interval:
X = zeros(1,K);
x = linspace(0,1,N);                                        %Interval
L = (abs(max(x))+abs(min(x)))/2;                            %Length
Fs = 8000;                                                 %SampleFreq
sek = N/Fs;                                                 %Reclength

%Functions:
%func = sin(2*pi*x);                                        %Sinus
%func = abs(sin(2*pi*x));                                   %Abs(Sinus)
%func = ones(1,N); func(1:floor(N/2)) = -ones(1,floor(N/2));%Wave puls
%func = x-0.5;                                               %Sawtooth
%func = x.^2;                                               %Square
%func = sin(2*pi*x)+x;                                      %x+sine
%func = exp(x);                                             %Exponential
%func = 2*pi*x.*sin(4*pi*x);                                %x*sin
%func = x.^x;                                               %Wierd
%func = zeros(1,N); func(floor(N/2)) = 1;                   %Dirac Delta
%func = 1/rand(1,1)*rand(1,N);                              %Random
%func = sort( 1/rand(1,1)*rand(1,N));                       %Sorted random
func = (record(sek,Fs))';                                  %Record

%Calculating average (two ways; a and b)
a = 0;
for n = 1:1:N-1
    a =a+((func(n+1)+func(n))/2*(x(n+1)-x(n)));
end
b = sum(func)/N;

%Calculating coeffs:
for k = 1:1:K;
    for n = 0:1:N-1
        X(k) = (X(k) + func(n+1)*exp(-1i*2*pi*k*(n+1)/N));
    end
end

%Calculating estimated function:
f=0;
g=0;
for k = 1:1:K
    f = f + 2/N*real(X(k)*exp(1i*k*pi*x/L));
    g = g + 2/N*imag(X(k)*exp(1i*k*pi*x/L));
end

% More periodes:
y = linspace(0,4,4*N);
Ly = (abs(max(y))+abs(min(y)))/2;
ff=0;
gg=0;
for k = 1:1:K
    ff = ff + 2/(4*N)*real(X(k)*exp(4*1i*k*pi*y/Ly));
    gg = gg + 2/(4*N)*imag(X(k)*exp(4*1i*k*pi*y/Ly));
end
X = [a,X];

%Plotting relevant data:
figure(1)
subplot(3,3,1),plot(x,func), title('Func(orig)'),ylabel('Amplitude'),xlabel('time'),axis([min(x)-0.5*max(x),max(x)+0.5*max(x),min(func)-0.5*max(func),max(func)+0.5*max(func)]);
subplot(3,3,2),plot(x,(a+f)), title('Func(esti)'),ylabel('Amplitude'),xlabel('time'),axis([min(x)-0.5*max(x),max(x)+0.5*max(x),min(a+f)-0.5*max(a+f),max(a+f)+0.5*max(a+f)]);
subplot(3,3,3),plot(x,func-(a+f)), title('Error'),ylabel('Amplitude'),xlabel('time'),axis([min(x)-0.5*max(x),max(x)+0.5*max(x),min(func-(a+f))-10*max(func-(a+f)),max(func-(a+f))+10*max(func-(a+f))]);
subplot(3,3,4),bar(real(X)),title('Real(coef)'),ylabel('Magnitude'),xlabel('frequency'), axis([-k/10,K*1.1,0,max(real(X))*1.1])
subplot(3,3,5),bar(imag(X)),title('Imag(coef)'),ylabel('Magnitude'),xlabel('frequency'), axis([-k/10,K*1.1,0,max(imag(X))*1.1])
subplot(3,3,6),plot(x,a+g),title('Imag(func)'),ylabel('Amplitude'),xlabel('time'),axis([min(x)-0.5*max(x),max(x)+0.5*max(x),min(a+g)-0.5*max(a+g),max(a+g)+0.5*max(a+g)]);
subplot(3,3,7),plot(y,a+ff),title('Real(esti)'),ylabel('Amplitude'),xlabel('time'),axis([min(y)-0.5*max(y),max(y)+0.5*max(y),min(a+ff)-0.5*max(a+ff),max(a+ff)+0.5*max(a+ff)]);
subplot(3,3,8),plot(y,a+gg),title('Imag(esti)'),ylabel('Amplitude'),xlabel('time'),axis([min(y)-0.5*max(y),max(y)+0.5*max(y),min(a+gg)-0.5*max(a+gg),max(a+gg)+0.5*max(a+gg)]);
subplot(3,3,9),plot(x,f-g),title('Func(Funn)'),ylabel('Amplitude'),xlabel('time'),axis([min(x)-0.5*max(x),max(x)+0.5*max(x),min(f-g)-0.5*max(f-g),max(f-g)+0.5*max(f-g)]);

%Write out error:
Error = norm(func-(a+f))
 
%Play recording:
playit((a+f),Fs);
pause(0.1);
playit(func,Fs);