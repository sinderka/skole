function A = K2(n)
A = spdiags(ones(n-1,1)*[-1 1],[0 1],n-1,n);
end