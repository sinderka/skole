%function [z] = intmethtest(zn,H,Hn,len_T,h,M)
function [y] = intmethtest(zn,A,v,len_T,h,M,Hn)
%%% Løser ligningen z'(t)-Hm*z = zn*v
%%%%%%%Skrive mer hva denne funksjonen gjør!!!!!!!!!!!
y = zeros(length(A),len_T);     % Skal bort
%z = zeros(M,len_T);
diff = 1;              %Skal bort

while (diff > 0.001)   %Skal bort
z = zeros(M,len_T);
   [V,H] = Arnoldi(A,v,M); %%% Skal bort
%[V,H,Hm] = krylov(A,v,M);
for j = 2:1:len_T
    z(:,j) = z(:,j-1) + H(1:M,1:M)*h*z(:,j-1);
    z(1,j) = z(1,j) + h*Hn*zn(j-1);
end


        zn = z(M,:);    %Skal bort
Hn = H(M+1,M);          % Skal bort
%Hn = Hm;                % Skal bort
nextstep = V(:,1:M)*z;  % Skal bort
y = y + nextstep;       % Skal bort

v = V(:,end);           % Skal bort
diff = max(max(nextstep));        % Vil ha den relativ!
end                         % Skal bort
end