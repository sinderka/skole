function coefff =  coeff(tri,p,i,k)
n = size(tri);
% k er trekant(1:n(1))
% i er hjørne (1:n(2))
% tri og p er punkter
temp = zeros(n(2),1);
temp(i) = 1;


coefff =     [p(tri(k,1),1),p(tri(k,1),2),1;
             p(tri(k,2),1),p(tri(k,2),2),1;
             p(tri(k,3),1),p(tri(k,3),2),1]\temp;
         
end
