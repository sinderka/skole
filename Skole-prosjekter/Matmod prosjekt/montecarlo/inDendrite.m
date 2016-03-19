function [ out ] = inDendrite( point )
%INDENDRITE Summary of this function goes here
%   Detailed explanation goes here
    
    radius = 0.22e-6;
    bodyHeight = 1e-6;

    cleftHeight = 15e-9;
    height = cleftHeight + bodyHeight + radius;

    % Check whether we exist in the main dendrite body (sylinder)
    inMainBody = point(2) >= height - bodyHeight && (point(1)-radius/2)^2 + (point(3) - radius/2)^2 <= radius^2;
    
    inTopHalfThingy =  (point(1)-radius/2)^2 + (height - point(2) - radius/2 - bodyHeight)^2 + (point(3) - radius/2)^2 <= radius^2;
    
    out = inMainBody || inTopHalfThingy;
end

