%%% For making input data,for the function "alf".
%%% Plotting out-data from the function "alf".


a = 0;

if a == 0
%l = [4,2,1,1,1,2,4,8]; lambda = l; m = l;my = 2; % length of the different links
%l = [1;1;2;1;4];lambda = l; m = l;my = 2;
%l = [1;1;5;1;1];lambda = l; m = l;my = 2;
%l = [0.5;0.56;0.5];lambda = l; m = l;my = 1;
%l = [1,2,3,4]; m = [5,6,7,8]; lambda = [9,10,11,12]; my = 13;
l = ones(1,20); m = l; lambda = l; my = 2; len = length(l);
%l = 0.05*ones(1,40); lambda = l;density = 1; m =density*l;my = 2;
%l = 10*rand(1,5); lambda = l;density = 1; m =rand(1,5);my = 2; len = length(l); m(3) = 100;
x = rand(len+1,1);
y = rand(len+1,1);
x(1) = 0;y(1) =0;       % Startpoints
x(end) = 1;y(end)=0;  % Endpoints
elseif (a == 1)
%%%Test1
l = [0.4,0.3,0.25,0.2,0.4]; lambda = l; m = l;my = 2;len = length(l);

x = rand(len+1,1);
y = rand(len+1,1);
x(1) = 0;y(1) =0;       % Startpoints
x(end) = 1;y(end)=-0.3;  % Endpoints
%%%Test2
elseif (a == 2)
l = [0.5,0.5,2,0.4,0.4]; lambda = l; m = l;my = 2;len = length(l);

x = rand(len+1,1);
y = rand(len+1,1);
x(1) = 0;y(1) =0;       % Startpoints
x(end) = 0;y(end)=-1;  % Endpoints
%%%Test3
elseif (a == 3)
l = [1,1]; lambda = l; m = l;my = 2;len = length(l);

x = rand(len+1,1);
y = rand(len+1,1);
x(1) = 0;y(1) =0;       % Startpoints
x(end) = 2;y(end)=0;  % Endpoints
%%%Test4
elseif (a == 4)
l = [1,1]; lambda = l; m = l;my = 2;len = length(l);

x = rand(len+1,1);
y = rand(len+1,1);
x(1) = 0;y(1) =0;       % Startpoints
x(end) = 0;y(end)=-2;  % Endpoints
end
len = length(l);
g = 9.81;


%x = rand(len+1,1);
%y = rand(len+1,1);
%x(1) = 0;y(1) =0;       % Startpoints
%x(end) = 1;y(end)=0;  % Endpoints

x = [x;y];
tic;
[x] = alf(l,x,lambda,my,g,m);
toc


%%% Plotting
x = reshape(x,length(x)/2,2);
plot(x(:,1),x(:,2),'*-')
%axis([-1 2 1 -10])
%axis([-0.5 0.5 -13 1])
