function [ out ] = findPhi(p,i)

% p er punkter, en matrise med [xvector,yvector]
% i er hvilket hjørne som phi skal finnes til.
n = size(p);

temp = zeros(n(1),1);

temp(i) = 1;


out = [p,ones(n(1),1)]\temp;

end

