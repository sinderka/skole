close all
clear all

kappa = 0.1;

k1 = 0.2;
k2 = 0.2;
k3 = 0.5;
k4 = 0.5;
k5 = 1;


T = 2;
dt = 0.1;

%h = 0.25;
h = 1;

U = zeros(7,50);
U(1,1) = 1;
U(2:7,1) = 0;


N = 5;


M = diag(ones(1,N - 2)*2 * h / 3) + diag(ones(1,N - 3) * h / 6, -1) + diag(ones(1,N - 3) * h / 6, 1);
Mhat = zeros(3 * (N - 2));
Mhat(N - 1 : 2*(N-2), N - 1 : 2*(N-2)) = eye(N-2);
Mhat(2*(N-2) + 1:3 * (N-2), 2*(N-2) + 1:3 * (N-2)) = eye(N-2);
Mhat(1 : N - 2 , 1 : N - 2) = M;

K = diag(ones(1,N - 2)* 2 / h ) + diag(ones(1,N - 3) * -1 / h , -1) + diag(ones(1,N - 3) * -1 / h, 1);
Khat = zeros(3 * (N - 2));
Khat(1 : N - 2 , 1 : N - 2) = M;

PaT = zeros(1, T/dt);
PaT(1) = 1;
QaPT = zeros(N);
QaPT(1,1) = 1 / kappa * PaT(1);

PbR = zeros(1, T/dt);
PbR(end) = 1;
QbPR = zeros(N);
QbPR(N,N) = 1 / kappa * PbR(end);

da = zeros(N,1);
da(1) = 1 / kappa;
daHat = [da ; 0 ; 0];


db = zeros(N,1);
db(N) = 1 / kappa;
dbHat = [db ; 0 ; 1];


ea = zeros(N,1);
eahat = [ea ; 1 ; 0];




% Xt = [U; Pa; Pb];
% Mhat = [M, zeros(N,2); zeros(2,N), [1, 0; 0, 1]];
% Khat = [K, zeros(N,2); zeros(2,N+2)];
% Qahat = [Qa, zeros(N,2); zeros(2,N), [U(1), 0; 0, 0]];
% Qbhat = [Qb, zeros(N,2); zeros(2,N), [0, 0; 0, U(N)]];
% dahat = [da; 0; 0];
% dbhat = [db; 0; 1];
% eahat = [zeros(N,1); 1; 0];

options = odeset('Mass', Mhat);

