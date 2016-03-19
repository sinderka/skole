function coefff =  coeff3(tri,p,i,k)
%n = size(tri);
% k er trekant(1:n(1))
% i er hjørne (1:n(2))
% tri og p er punkter
%temp = eye(4);


%coefff =     [p(tri(k,1),1),p(tri(k,1),2),p(tri(k,1),3),1;
%              p(tri(k,2),1),p(tri(k,2),2),p(tri(k,2),3),1;
%              p(tri(k,3),1),p(tri(k,3),2),p(tri(k,3),3),1;
%              p(tri(k,4),1),p(tri(k,4),2),p(tri(k,4),3),1]\temp;
%                 
%
temp = ones(4,1);
coefff = [p(tri(k,:),:),temp]\eye(4);         
%coefff = coefff(i,:)';     Eller
coefff=coefff(:,i);
end
