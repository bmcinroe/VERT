function [nPlane, nExplained, max_tangent_dim, cov_matrix] = construct_vielbein(D_i, neighborhood_i, thresh)

nPlane = [];
nExplained = [];

[coeff,score,latent,tsquared,explained,mu] = pca([D_i ; neighborhood_i]); % PCA of point and its neighbors

% Generate nx1 structs for VB's and % variances explained up to tangent_dim
tangent_dim = 1;

while explained(tangent_dim) > thresh
    nPlane = [nPlane, coeff(:,tangent_dim)];
    nExplained = [nExplained; sum(explained(1:tangent_dim))];
    
    tangent_dim = tangent_dim + 1;
    if tangent_dim == length(explained)
        break
    end
end

max_tangent_dim = tangent_dim - 1;

cov_matrix = cov([D_i ; neighborhood_i]);

end