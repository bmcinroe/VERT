function [nIndex, neighborhoods] = find_neighborhoods(D, eps, partitions, window)

neighborhoods = cell(length(partitions), 1);
nIndex = cell(length(partitions), 1);

% Local Search

for i = window(1):window(2)
    euclidnorms = vecnorm(D(partitions{i}, :) - D(i,:), 2, 2);
    neighbors = find(euclidnorms < eps);
    
    nIndex{i} = partitions{i}(neighbors); % indicies of neighbors
    neighborhoods{i} = D(partitions{i}(neighbors), :);
end


end

