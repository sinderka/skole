function k = getTriangle3(tri,p1,p2,p3)
for i = 1:max(size(tri))
    for j = 1:min(size(tri))
        for l = 1:min(size(tri))
            for m = 1:min(size(tri))
                if tri(i,j) ==p1 && tri(i,l) == p2 && tri(i,m) == p3
                    k(1) = i;
                    k(2) = j;
                    return;
                end
            end
        end
    end
end

error('Triangle not found!')
end