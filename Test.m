%close all;
k = 8;          % Ikke skift før det fungerer

f = @(x)1;
%f = @(x) sin(x)*x-x^x;
A = diag(1:k);
%A =rand(k,k);
%u = rand(k,1);
u  = (1:k)';
%A = gallery('tridiag',k,1,4,1);
len_v = k;
len_T = 10000;
T = linspace(0,2,len_T);


v = A*u;
M = 2;      % Ikke skift før det fungerer

%y = zeros(len_v,len_T);
%%%%%%%%%%%%%%<Første runde>%%%%%%%%%%%
zn = zeros(M,len_T);
[Vn,Hn] = Arnoldi(A,v,M);

for j = 2:1:len_T
    h = T(j)-T(j-1);
    zn(:,j) = zn(:,j-1)+h*Hn(1:M,1:M)*zn(:,j-1);
    zn(1,j) = zn(1,j) + h*f(T(j-1))*norm(v,2);
end

y = Vn(:,1:M)*zn; %Arg, føler at jeg gjør bare er feil hele tiden!

%y =y+ yn;           % Burde det være noen ristriksjoner på dimensjoner eller noe slikt?
%%%%%%%%%%%%%%%%%%%%%%%%%%plotteting
figure(1);
sol = zeros(length(v),len_T);
for i = 1:1:len_T
    sol(:,i) = (expm(A*T(i))-eye(k))*u;
end
plot(T,sol,'b')
hold on
plot(T,y,'g')
hold off
for i = 1:k
    Err = norm(sol(i,:)-y(i,:) ,Inf)/max(abs(sol(i,:)));
    fprintf('Relative error for graph %d:\t ||y-y_n||_Inf = %f\n',i, Err )
end
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%plotteting
%%%%%%%%%%%%%</Første runde>%%%%%%%%%%%


% Men første runde fungerer fint!



vn = Vn(:,end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%<Andre runde>%%%%%%%%%%%
zm = zeros(M,len_T);
[Vm,Hm] = Arnoldi(A,vn,M);
%[Vm,Hm] = Arnoldi(A,vn,M);
vm = Vm(:,end);

for j = 2:1:len_T
    h = T(j)-T(j-1);
    %err = Hn(M+1,M)*zn(M,j-1)*v; %%% Virker ikke som den rette(feilen økte)
    err = Hn(M+1,M)*zn(M,j-1)*vn; %%%DENNE
    %err = Hn(M+1,M)*zn(M,j-1)*vm; %%%Ikke testet
    
    %err = Hn(M+1,M)*zm(M,j-1)*v; %%%Feilen forandret seg ikke
    %err = Hn(M+1,M)*zm(M,j-1)*vn; %%%Feilen forandret seg ikke
    %err = Hn(M+1,M)*zm(M,j-1)*vm; %%%Ikke testet
    
    %err = Hm(M+1,M)*zm(M,j-1)*v; %%%Feilen forandret seg ikke
    %err = Hm(M+1,M)*zm(M,j-1)*vn; %%%ELLER DENNE
    %err = Hm(M+1,M)*zm(M,j-1)*vm; %%%Ikke testet
    
    %err = Hm(M+1,M)*zn(M,j-1)*v; %%%Feilen forandret seg ikke
    %err = Hm(M+1,M)*zn(M,j-1)*vn;
    %err = Hm(M+1,M)*zn(M,j-1)*vm; %%%Ikke testet
    
    zm(:,j) = zm(:,j-1) + h*Vm(:,1:M)'*(A*Vm(:,1:M)*zm(:,j-1)+err); %Den

    %zm(:,j) = zm(:,j-1)+h*(Vm(:,1:M)'*A*Vm(:,1:M)*zm(:,j)+Vm(:,1:M)'*zn(M,j-1)*v*Hn(M+1,M));
end



y = y+Vm(:,1:M)*zm; % Skal Vn fnuttes
%y =y+ ym;
%y = ym;

%%%%%%%%%%%%%%%%%%%%%%%%plotteting
sol = zeros(length(v),len_T);
for i = 1:1:len_T
    sol(:,i) = (expm(A*T(i))-eye(k))*u;
end
figure(2);
plot(T,sol,'b')
hold on
plot(T,y,'g')
hold off


for i = 1:k
    Err = norm(sol(i,:)-y(i,:) ,Inf)/max(abs(sol(i,:)));
    fprintf('Relative error for graph %d:\t ||y-y_n||_Inf = %f\n',i, Err )
end
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%plotteting
%%%%%%%%%%%%%%</Andre runde>%%%%%%%%%%%



v = Vm(:,end);


%%%%%%%%%%%%%%<Tredje runde>%%%%%%%%%%%
[Vk,Hk] = Arnoldi(A,v,M);



zk = zeros(M,len_T);
%temp = zeros(len_v,1);

for j = 2:1:len_T
    %h = T(j)-T(j-1);
    %zk(:,j) = zk(:,j-1)+h*(Vk(:,1:M)'*A*Vk(:,1:M)*zk(:,j)+Vk(:,1:M)'*zm(M,j-1)*v);
    
    
    %%%
    h = T(j)-T(j-1);
    %err = Hn(M+1,M)*zn(M,j-1)*v; %%% Virker ikke som den rette(feilen økte)
    %err = Hn(M+1,M)*zn(M,j-1)*vn; %%%DENNE
    %err = Hn(M+1,M)*zn(M,j-1)*vm; %%%Ikke testet
    
    %err = Hn(M+1,M)*zm(M,j-1)*v; %%%Feilen forandret seg ikke
    %err = Hn(M+1,M)*zm(M,j-1)*vn; %%%Feilen forandret seg ikke
    %err = Hn(M+1,M)*zm(M,j-1)*vm; %%%Ikke testet
    
    %err = Hm(M+1,M)*zm(M,j-1)*v; %%%Feilen forandret seg ikke
    %err = Hm(M+1,M)*zm(M,j-1)*vn; %%%ELLER DENNE
    %err = Hm(M+1,M)*zm(M,j-1)*vm; %%%Ikke testet
    
    %err = Hm(M+1,M)*zn(M,j-1)*v; %%%Feilen forandret seg ikke
    %err = Hm(M+1,M)*zn(M,j-1)*vn;
    err = Hm(M+1,M)*zm(M,j-1)*vm; %%%Ikke testet
    
    zk(:,j) = zk(:,j-1) + h*Vk(:,1:M)'*(A*Vk(:,1:M)*zk(:,j-1)+err); %Den
    %%%
end



y =y+Vk(:,1:M)*zk; % Skal Vn fnuttes
%y =y+ ym;
%y = yk;
%%%%%%%%%%%%%%</Tredje runde>%%%%%%%%%%%

figure(3);
plot(T,sol,'b')
hold on
plot(T,y,'g')
hold off


for i = 1:k
    Err = norm(sol(i,:)-y(i,:) ,Inf)/max(abs(sol(i,:)));
    fprintf('Relative error for graph %d:\t ||y-y_n||_Inf = %f\n',i, Err )
end
fprintf('\n')




