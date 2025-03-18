function partitions = partition_samples(D, rho, window)

% window: a 2 X 1 array containing the initial and final index for analysis

partitions = cell(length(D), 1); % cell array of partitions for each sample

[coeff,score,latent,tsquared,explained,mu] = pca(D);

theta = D * coeff(:,1);

for i = window(1):window(2)
    partitions{i} = find(theta < theta(i) + rho & theta > theta(i) - rho);
end

end

