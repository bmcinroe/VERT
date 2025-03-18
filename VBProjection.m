function VBProjections = VBProjection(Displacements, VB)

VB_dim = size(VB, 2);

VBProjections = zeros(size(Displacements, 1), VB_dim);

for i = 1:size(Displacements, 1)
    
    VBProjections_i = zeros(1, VB_dim);
    
    for j = 1:VB_dim
        
        VBProjections_i(j) = dot(Displacements(i,:), VB(:,j)) ./ dot(VB(:,j), VB(:,j));
        
    end
    
    VBProjections(i,:) = VBProjections_i;

end

end

