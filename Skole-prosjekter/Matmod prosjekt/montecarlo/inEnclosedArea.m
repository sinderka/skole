function [ outBool ] = inEnclosedArea( diffs )
%INENCLOSEDAREA Summary of this function goes here
%   This function defines our enclosed area. No particle can move out of
%   here

    radius = 0.22e-6;
    bodyHeight = 1e-6;

    cleftHeight = 15e-9;
    height = cleftHeight;% + bodyHeight + radius;


inMainArea = ~(diffs(1) < 0  || diffs(1) > radius || diffs(2) < 0 || diffs(2) > height || diffs(3) < 0 || diffs(3) > radius);

%inLeftCircle = sqrt((diffs(1)/radius)^2 + (diffs(2)/height)^2) < 1;
outBool = inMainArea;% && ~inLeftCircle;

end

