close all;
k = 25;

%k antall interne noder (i hver rettning)
%k+2 antall noder totalt (i hver rettning)
% steglengde  er da 1/k(+1)



%f1 = @(x) 5; %f1 = @(x)5*x;
%f2 = @(x)3;%f2 = @(x)12*x^2;
%f = @(x) sin(x)*x-x^x;
f1 = @(t) 1/(t+1)^2;
f2 = @(t) t/(t+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%Få hjelp til å se hva som bruker mest tid!%%%
%A =rand(k,k);
hs = 1/(k+1);
%A = hs^2*gallery('poisson', k);
A = gallery('tridiag',k,-1,2,-1);
%u = rand(k,1);
%u1 = rand(k,1);
%u2 = rand(k,1);
funcU1 = @(x) x.*(x-1);
funcU2 = @(x,y) -2;
Y = linspace(0,1,k+2)'; %X = linspace(0,1,k+2)';
%u1 = zeros(k,1); u2 = zeros(k,1);       % Burde jeg ta med kanten her?
%for i = 1:k
   u1 = funcU1(Y(2:k+1));
   u2 = funcU2(Y(2:k+1))*ones(k,1);
%end




%Y = linspace(0,1,k); X = linspace(0,1,k);

%A = gallery('tridiag',k,1,4,1);
%u1 = 1:1:k; u2 = -1:-1:-k; u1 = u1'; u2 = u2';
len_T = 100; T = linspace(0,2,len_T);
h = T(2)-T(1);

%v = A*u;
v1 = A*u1;
v2 = A*u2;
%v1 = u1;
%v2 = u2;

disk =h/hs^2; %burde være liten (<1/2)





M = 1;  %Restarte variable
%[Va,Ha]=Arnoldi(A,v1,M)
%[Vk,Hk,Hm]=krylov(A,v1,M)

zn1 = zeros(1,len_T);
zn2 = zeros(1,len_T);
for i=1:1:len_T
zn1(i) = f1(T(i));
zn2(i) = f2(T(i));
end



%%% Forskjellige implementasjoner
%[y] = resArnoldi(A,v1,f1,M,T);
%y = yn;

%[y] = intmethtest(zn1,A,v1,len_T,h,M, norm(v1,2));
%[y] = KPM_IVP_res(A,v1,M,zn1,h,len_T);



%%%%% Ting som burde fungere
y = 0;
diff =1;
Hn1 = norm(v1,2); Hn2 = norm(v2,2);
while diff > 0.001 %%%Finne en bedre måte å rengne ut feilen på?
%[V,H] = Arnoldi(A,v1,M);
[V1,H1,Hm1] = krylov(A,v1,M);
[V2,H2,Hm2] = krylov(A,v2,M);

[zn1] = intmeth(zn1(end,:),H1,Hn1,len_T,h,M,1);
[zn2] = intmeth(zn2(end,:),H2,Hn2,len_T,h,M,1);
%[zn1] = intmeth(zn1(end,:),H,Hn,len_T,h,M,2);
%[zn1] = intmeth(zn1(end,:),H,Hn,len_T,h,M,3);
%Hn = H(M+1,M);
Hn1 = Hm1; Hn2 = Hm2;
ns = V1(:,1:M)*zn1+V2(:,1:M)*zn2;
y = y + ns;
v1 = V1(:,end); v2 = V2(:,end);
diff = max(max(ns));
end


solution = @(t,x) t*x.*(x-1)/(t+1);


% ta med nullene på kanten:
sol = zeros((k+2),len_T);
for i = 1:len_T
    %for j = 1:k+2
            sol(:,i) = solution(T(i),Y);
    %end
end

% %Ikke ta med nullene på kanten
% sol = zeros(k^2,len_T);
% for i = 1:len_T
%     for j = 1:k
%             sol(1+(j-1)*(k):j*(k),i) = solution(T(i),X(j+1),Y(2:end-1));
%     end
% end

        
%%% legge til masse nuller på y :(
% ynew = zeros(k+2,len_T);
% ynew(2:end-1,:) = y;
% y = ynew;



%sol = zeros(k^2,len_T);
%for i = 2:1:len_T %Burde lage en metode som regner ut "korrekt" løsning med differensementoder
    %sol(:,i) = (expm(A*T(i))-eye(k))*(u1*f1(1)+u2*f2(1));
%    sol(:,i) = (expm(A*T(i))-eye(k))*(u1*f1(1));
    %sol(:,i)= sol(:,i-1) + h*A*sol(:,i-1)+ h*f1(T(i-1))*v1; %% Hvorfor er ikke dette rett? Ingenting fungerer jo lenger!
%end

%%% burde fikse slik at jeg kan se hva som skjer!
%for i = 1:10:len_T
%ytemp = reshape(y(:,i),[k+2,k+2]);
figure(1)
surf(y(2:end-1,:))
%surf(y)
%surf(X,Y,ytemp)
%drawnow
%mesh(X,Y,y)
%hold on
%plot(T,sol,'b')
%soltemp = reshape(sol(:,i),[k+2,k+2]);
figure(2)
surf(sol)
%drawnow
%end
%figure(3)
%scatter(X,Y,sol)
%figure(3)
%plot(T,sol-y,'b')
%figure(4)
%plot(T,y,'g')


Err = max(max(abs(sol-y)))/max(max(abs(sol)));
fprintf('Relative error:%f\n', Err )
%for i = 1:k
%    Err = norm(sol(i,:)-y(i,:) ,Inf)/max(abs(sol(i,:)));
%    fprintf('Relative error for graph %d:\t ||y-y_n||_Inf = %f\n',i, Err )
%end
