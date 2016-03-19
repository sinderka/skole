




%%%%%%%%%%%%%Start på å lage input data%%%%%%%%%%%%%%%%%%%%%
m = 1000;

f1 = @(t) 1/(t+1).^2;                        %1
f2 = @(t) t-t+1;                               %2


hs = 1/(m+1);
A = hs^2*gallery('poisson', m);


Y = linspace(0,1,m+2)'; X = Y;

g1 = @(x,y) x*(x-1)*y.*(y-1);               %1
g2 = @(x,y) x*(x-1)*sin(pi*y);              %2


v1 = zeros(m^2,1); v2 = zeros(m^2,1);
for i = 1:m
    v1(1+(i-1)*m:i*m) = g1(X(i+1),Y(2:m+1));    %1
    v2(1+(i-1)*m:i*m) = g2(X(i+1),Y(2:m+1));    %2
end

k = 100; T = linspace(0,1,k); ht = T(2)-T(1);


n = 2;                              % restart variable
intmet = 2;  %Sier noe om intgrasjonsmetode, sjekk "intmeth.m" for mer informasjon

%%% Time dependancy:
F1 = zeros(1,k) ; F2 = zeros(1,k);
for i = 1:k
    F1(i) = f1(T(i));
    F2(i) = f2(T(i))';
end
%%%%%%%%%%%%%Slutt på å lage input data%%%%%%%%%%%%%%%%%%%%%
%%% Faktisk metode %%%
U1 = projmeth(F1,A,v1,k,ht,n,intmet);
U2 = projmeth(F2,A,v2,k,ht,n,intmet);
U = U1+U2;
%%% Faktisk metode %%%


%%% Sammenligne stuff %%%
compare = 1;
if compare
    %%%% Solution %%%%
solution1 = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1);             %1
solution2 = @ (t,x,y) t*x*(x-1)*sin(pi*y);                  %2
solution = @(t,x,y) solution1(t,x,y) + solution2(t,x,y);    %0
sol = zeros(m^2,k);
%%% Test funksjon
for i = 1:k
    for j = 1:m
        sol(1+(j-1)*(m):j*(m),i) = solution(T(i),X(j+1),Y(2:end-1));
    end
end
%%%% Solution %%%%
%%% Ta med nuller på kanten:


solk = zeros((m+2)^2,k);
for i = 1:k
    for j = 1:m
        solk(1+(j-1)*(m+2):j*(m+2),i) = solution(T(i),X(j),Y);
    end
end
%Ikke ta med nullene på kanten:

%%% Konstant funksjon
% for i = 2:1:len_T %Burde lage en metode som regner ut "korrekt" løsning med differensementoder
%    sol(:,i) = (expm(A*T(i))-eye(k^2))*(v1*f1(1)+v2*f2(1)); % Om f1==const
%    %sol(:,i) = (expm(A*T(i))-eye(k))*(v1*f1(1)); % Om f==const
% end




Numsol = zeros(m^2,k);
%Numsol(:,1) = init; %Om startbetingelser
for i = 2:1:k %Burde lage en metode som regner ut "korrekt" løsning med differensementoder
    %sol(:,i) = (expm(A*T(i))-eye(k^2))*(v1*f1(1)+v2*f2(1)); % Om f1==const
    %sol(:,i) = (expm(A*T(i))-eye(k))*(v1*f1(1)); % Om f==const
    Numsol(:,i)= Numsol(:,i-1) + ht*A*Numsol(:,i-1)+ ht*f1(T(i-1))*v1+ht*f2(T(i-1))*v2; %% Om f ~= const
    %Numsol(:,i)= Numsol(:,i-1) + h*A*Numsol(:,i-1)+ h*f1(T(i-1))*v1+h*f2(T(i-1))*v2; %% Om f ~= const
end

%%% Endre dimensjon på y og legge til nuller på randen
% film = 0;
% if film
% ynew = zeros((m+2),m+2,k);
% ry = reshape(U,[m,m,k]);
% ynew(2:end-1,2:end-1,:) = ry;
% solnew = reshape(solk,[m+2,m+2,k]);
% 
% for i = 1:1:k
%     figure(5)
%     surf(X,Y,ynew(:,:,i))
%     drawnow
%     figure(6)
%     surf(X,Y,ynew(:,:,i)-solnew(:,:,i))
%     drawnow
%     
%     
% end
% end


% figure(1)
% plot(T,U,'g')
% figure(2)
% plot(T,sol,'b')
% figure(3)
% plot(T,Numsol,'r');
%figure(2)
%surf(sol)
%figure(3)
%plot(T,sol,'b')
%figure(4)
%plot(T,y,'g')


%Err = max(max(abs(sol-y)))/max(max(abs(sol)));
%fprintf('Relative error:%f\n', Err )
Err = zeros(1,m);
for i = 1:m
    Err(i) = norm(sol(i,:)-U(i,:) ,Inf)/max(abs(sol(i,:)));
    %fprintf('Relative error for graph %d:\t ||y-y_n||_Inf = %f\n',i, Err(i) )
end
fprintf('Max err: %f \n', max(Err));
end