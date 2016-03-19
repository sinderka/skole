% Inspiration found at: http://people.sc.fsu.edu/~jburkardt/m_src/brownian_motion_simulation/brownian_motion_simulation.html
close all
clear all

M = 3;              % dimensions:
N = 100;            % timesteps
t = 2e-2;        % time at N't timestep
diff = 0.3e-12;  % diffusion coefficient
NN = 5000;          % number of particles

radius = 0.22e-6;
height = 15e-9;
numOfReceptors = 1000*pi*(radius*1e6)^2;

dt = t/(N-1);
ts = 0:dt:t;

coords = zeros(3*NN,N);
coords(1:3:end,1) = 0.5*radius;        % Skewing to middle of x
coords(3:3:end,1) = 0.5*radius;        % Skewing to middle of z

cleanup = false;
% store statuses to all neurotransmitters:
status = zeros(NN,N);

% Statuses:
% 0 = free, open thingy
% 1 = stuck at top edge
% 3 = stuck at GLIA cells
% 4 = free, deactivated
% 5 = collected at axon again

%%

for j = 2:N
    for i = 1:NN
        if ~cleanup && sum(status(:,j) == 1) >= 0.8*numOfReceptors
            % Start cleanup process
            cleanup = true;
            TIME = ts(j);
        end
        
        point = coords(3*i-2:3*i, j-1);
        
        % Generating random step
        length = sqrt ( diff * dt ) * randn ( 1 );
        direction = randn(M,1);
        step = length*direction;
        newPoint = point + step; 
        while ~inEnclosedArea(newPoint)
            % if destination is outside area, generate new step
            length = sqrt ( diff * dt ) * randn ( 1 );
            direction = randn(M,1);
            step = length*direction;
            newPoint = point + step;
        end
        
        z = rand(1);
        switch status(i,j)
            case 0  % free activated
                if cleanup == false
                    % Checking the density of receptors around this area:  
                    occupied = sum(status(:,j) == 1);
                    percentAvailable = (numOfReceptors-occupied)/numOfReceptors;
                    
                    percentToWall = (height - newPoint(2))/height;
                    if (percentToWall < 0.1)*(1-percentToWall)*percentAvailable >= z*20
                        newPoint(2) = height;      %stuck at wall
                        status(i,j:end) = 1;
                    end
                else % we want to clean up with glia cells:
                    mypercentLeft = (radius-newPoint(1))/radius;
                    mypercentRight = newPoint(1)/radius;
                    
                    if (mypercentLeft < 0.1)*(1-mypercentLeft) >= z*20
                        newPoint(1) = radius; 
                        status(i,j:end) = 3;
                    elseif (mypercentRight < 0.1)*(1-mypercentRight) >= z*20
                        newPoint(1) = 0; 
                        status(i,j:end) = 3;
                    end
                    
                end
                
            case 1  % Stuck at dendrite wall
                if z >= 0.98 || cleanup      
                    status(i,j:end) = 0; 
                else
                    % We have a chance it is just stuck there:
                    newPoint = point;
                end
                
            case 3 % stuck at GLIA cell, just release:
                status(i,j:end) = 4;
                
            case 4 % free deacticvated                
                percentToWall = newPoint(2)/height;
                if (percentToWall < 0.1)*(1-percentToWall) >= z
                    newPoint(2) = 0;      %stuck at wall
                    status(i,j:end) = 5;
                end
                
            case 5 % stuck at axon:
                % Dont change anything
                newPoint = point;
                
        end
        % Applying changes to coordinate matrix
        coords(3*i-2:3*i, j) = newPoint;
         
    end 
    j
end
%%
figure
plot(ts,sum(status == 1), '*')
xlabel('Time');
ylabel('Bound receptors');
%title('number of stuck particles')

%%
fig = figure;
set(fig, 'Visible','off')
view(3)
for i = 1:N
    clf
    hold on
    view ( 3 )
    
    x = coords(1:3:end,i);
    y = coords(2:3:end,i);
    z = coords(3:3:end,i);
    
    % Plotting with different status:
    plot3 ( x(status(:,i) == 0), y(status(:,i) == 0), z(status(:,i) == 0), '.b')
    plot3 ( x(status(:,i) == 1), y(status(:,i) == 1), z(status(:,i) == 1), '*c')
    plot3 ( x(status(:,i) == 3), y(status(:,i) == 3), z(status(:,i) == 3), '*g')
    plot3 ( x(status(:,i) == 4), y(status(:,i) == 4), z(status(:,i) == 4), '.r')
    plot3 ( x(status(:,i) == 5), y(status(:,i) == 5), z(status(:,i) == 5), '*black')
    
    axis([0 radius 0 height 0 radius])
    grid on
    xlabel ( 'X' ), ylabel ( 'Y' ), zlabel ( 'Z' )
    
    title ( '3d Monte Carlo simulation' )
    saveas(fig, sprintf('../figures/time%03d.png',i));
    pause(0.1)
end
