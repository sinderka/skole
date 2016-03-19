function [y] = resArnoldi(A,v,f,M,T)


%close all;
%k = 5000;          % Ikke skift før det fungerer

%f = @(x)1;
%f = @(x) sin(x)*x-x^x;
%f = @(x) x^2+x*sin(cos(x+x^2)) + 3*x;
%A = diag(1:k);
%A =rand(k,k);
%u = rand(k,1);
%u  = (1:k)';
%A = gallery('tridiag',k,1,4,1);
%len_v = k;
%len_T = 1000;
%T = linspace(0,2,len_T);


%v = A*u;
%M = 4;



len_v = length(A);
len_T = length(T);

y = zeros(len_v,len_T);

diff = 1;
first = 1;

while first || (diff > 0.001 && M ~= len_v)
    
    
    z = zeros(M,len_T);
    [V,H] = Arnoldi(A,v,M);
    if first
        
                %%% Bruke matlabs ting for å fikse ting
                E_1 = zeros(M,1); E_1(1) = 1;
                [~,z] = ode45( @(t,y) H(1:M,1:M)*y+f(t)*norm(v,2)*E_1,T,zeros(M,1),['RelTol',1e-6]);
                z = z';
                
                %Bruker korrekt integrasjon, men da må f(t) == konst;
                %z = @(t) H(1:M,1:M)\(expm(H(1:M,1:M)*t)-eye(M))*f(1)*norm(v,2)*E_1;

                
         %for j = 2:1:len_T

                %Bruker korrekt integrasjon, men da må f(t) == konst;
         %       z(:,j) = H(1:M,1:M)\(expm(H(1:M,1:M)*T(j))-eye(M))*f(1)*norm(v,2)*E_1;
                
                
              %%%%Bruker frammover eller bakover differanser:
%             h = T(j)-T(j-1);
%             z(:,j) = z(:,j-1)+h*H(1:M,1:M)*z(:,j-1);
%             z(1,j) = z(1,j) + h*f(T(j-1))*norm(v,2);
% 
         %end
                    first = 0;
    else
        
            %%%% Andre ordens metode:
           %B =  gallery('tridiag',len_T,1,0,-1);
        
        %[T,z] = ode45( @(t,y) V(:,1:M)'*A*V(:,1:M)*y+V(:,1:M)'*err(t),T,[0,0],['RelTol',1e-3]);
        %tempmat = V(:,1:M)'*A*V(:,1:M);
        for j = 2:1:len_T
            %%%%Kanskje ikke rett fordi zn = z-zm eller noe slikt (Men jeg mener det er korrekt!)
            
            
            
            
            %Bruker korrekt integrasjon, men da må f(t) == konst; Fungerer
            %ikke enda ( Men er "lett" og fikse) 
            %z(:,j) = H(1:M,1:M)\(expm(H(1:M,1:M)*T(j))-eye(M))*zn(M,j)*norm(v,2)*E_1;
            %%%% 
            
            %%%% Fungerer ikke enda
            %z(:,j) = quadrature1D([
            %T(j),T(j-1)],[j,j-1],1,tempmat*z(:,j-1)+V(:,1:M)'*err);
            %Fungerer ikke for matriser
            
            %%% Ikke så veldig bra
            %[~,temp] = ode45( @(t,y) H(1:M,1:M)*y+zn(j)*norm(v,2)*E_1,[T(j),T(j-1)],zeros(M,1),['RelTol',1e-6]);
            %z(:,j) = temp(end,:)';
            
            %%%%Dårlig første ordens metode
             h = T(j)-T(j-1);
             z(:,j) = z(:,j-1) + h*(H(1:M,1:M)*z(:,j-1));
             z(1,j) = z(1,j) + h*Hn*zn(j-1);
        end
    end
    zn = z(M,:);
    Hn = H(M+1,M);
    nextstep = V(:,1:M)*z;
    y = y + nextstep;
    
    v = V(:,end);
    diff = max(max(nextstep));        % Vil ha den relativ!

end

% Neste steg: 
% 1. integrere på en mer fornuftig måte. Nope! 
% 2. Lagre konstanter/variabler på en mer kostnadseffektiv måte. JeeiXXX
% 3. Gjøre klart for paralellisering.
% 4. Forenkle err termen og slikt. Jeei, har forenklet litt XXX



% figure(1);
% sol = zeros(length(v),len_T);
% for i = 1:1:len_T
%    sol(:,i) = (expm(A*T(i))-eye(k))*u;
% end
% plot(T,sol,'b')
% hold on
% plot(T,y,'g')
% hold off
% for i = 1:k
%     Err = norm(sol(i,:)-y(i,:) ,Inf)/max(abs(sol(i,:)));
%     fprintf('Relative error for graph %d:\t ||y-y_n||_Inf = %i\n',i, Err )
 end