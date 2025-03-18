function [templateMetric, beta] = templateMetric(VB_Position_i, VB_Displacement_i, VB_Positions, VB_Displacements)

v_star = (VB_Displacement_i./norm(VB_Displacement_i)).';

% Center local vector field data
position_centered = VB_Positions - VB_Position_i;

displacement_centered = VB_Displacements - VB_Displacement_i;

beta = (pinv(position_centered.'*position_centered) * (position_centered.' * displacement_centered));

beta = beta.';
v_hat = beta * position_centered.';
v_hat = v_hat.';

[u,s,w] = svd(beta * beta.');

% Used for kernel metric - deprecated
%template_dim = 1
%template_subspace = w(:,(size(w,2) - template_dim + 1):size(w,2));
%templateMetric = subspace(v_star, template_subspace);

% inner product metric
if size(beta) == [0 0] 
    templateMetric = NaN; % return NaN if local Jacobian could not be estimated (e.g., no neighbors found)
else
    templateMetric = norm(beta * VB_Displacement_i.');
end

% FTLE Metric
% if size(alpha) == [0 0] 
%     templateMetric = NaN;
% else
%     templateMetric = 0.5 * log(max( eig( beta * beta.' ) ) );
% end

%template_subspace_proj = inv(template_subspace.' * template_subspace) * (template_subspace.' * v_star)

%templateMetric = acos(dot(v_star, template_subspace_proj)./(norm(v_star)*norm(template_subspace_proj)));

end

