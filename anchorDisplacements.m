function nDifferences = anchorDisplacements(D, nIndex, traj_I, traj_F, h)
% Computes finite differences for all samples within a neighborhood

% Conditional ensures that points at start and end of trajectory take
% forward and backward differences explicitly

% h: time step 

nDifferences = zeros(length(nIndex), size(D,2));

for i = 1:length(nIndex)
    
    idx = nIndex(i);
    
    if ismember(idx, traj_I)
        
        nDifferences(i,:) = (D(idx + 1,:) - D(idx,:))./h; % compute forward diff if IC
        
    elseif ismember(idx, traj_F)
        
        nDifferences(i,:) = (D(idx,:) - D(idx - 1,:))./h; % compute backward diff if FC
        
    else
        
        nDifferences(i,:) = (D(idx + 1,:) - D(idx - 1,:))./(2*h); % central diff otherwise
        
    end
    
        
end

% for i = 1:length(nIndex)
%     
%     idx = nIndex(i);
%     
%     if idx == window(1) || idx == 1
%         
%         nDifferences(i,:) = (D(idx + 1,:) - D(idx,:))./h;
%         
%     elseif idx == window(2) || idx == size(D,1)
%         
%         nDifferences(i,:) = (D(idx,:) - D(idx - 1,:))./h;
%         
%     else
%         
%         nDifferences(i,:) = (D(idx + 1,:) - D(idx - 1,:))./(2*h);
%         
%     end
%     
%         
% end


end