function [neighborhoods, nIndex, nSize] = kNN_neighbor_reduction(k, D, nIndex_euclid, neighborhoods_euclid, window)

neighborhoods = cell(length(nIndex_euclid), 1);
nIndex = cell(length(nIndex_euclid), 1);
nSize = zeros(length(nIndex_euclid), 1);

for i = window(1):window(2)
    Idx = knnsearch(neighborhoods_euclid{i}, D(i,:), 'K', k);
    
    neighborhoods{i} = neighborhoods_euclid{i}(Idx,:);
    nIndex{i} = nIndex_euclid{i}(Idx);
    nSize(i) = length(nIndex{i});
end

end

