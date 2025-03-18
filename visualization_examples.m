%% Example plotting methods
% grey color for visuals
grey = [224 224 224]./255; 
% GPC axes for data visualization
[coeff,score,latent,tsquared,explained,mu] = pca(D);

%% Visualize distance estimator over projected samples
% Estimator visualization against all samples
figure;
hold on;
title("Fiberwise Distance")
scatter3(s(:,:) * coeff(:,1), s(:,:) * coeff(:,2), s(:,:) * coeff(:,3),1, 'o', 'MarkerFaceColor', grey, 'MarkerEdgeColor', grey);
scatter3(s(window(1):window(2),:) * coeff(:,1), s(window(1):window(2),:) * coeff(:,2), s(window(1):window(2),:) * coeff(:,3), 20, [Vielbein.Metric], 'filled'); colormap jet; colorbar; 
hold off;

% Estimator Visualization against all samples (log version)
figure;
hold on; 
title("Fiberwise Distance (log)")
scatter3(s(:,:) * coeff(:,1), s(:,:) * coeff(:,2), s(:,:) * coeff(:,3),1, 'o', 'MarkerFaceColor', grey, 'MarkerEdgeColor', grey);
scatter3(s(window(1):window(2),:) * coeff(:,1), s(window(1):window(2),:) * coeff(:,2), s(window(1):window(2),:) * coeff(:,3), 20, log([Vielbein.Metric]), 'filled'); colormap jet; colorbar; 
hold off;

%% Visualize neighborhood
% Plot neighbors
figure;
hold on;
title("Neighborhood Visualization")
m=540; % Index of sample whose neighborhood will be plotted
scatter3(s(:,:) * coeff(:,1), s(:,:) * coeff(:,2), s(:,:) * coeff(:,3),1, 'o', 'MarkerFaceColor', grey, 'MarkerEdgeColor', grey);
scatter3(s([Vielbein(m).Neighbors],:) * coeff(:,1), s([Vielbein(m).Neighbors],:) * coeff(:,2), s([Vielbein(m).Neighbors],:) * coeff(:,3),80, 'r.');
hold off;

%% Visualize normalized estimator trace
% Plot normalized trace of distance estimator (blue) with representative subsystem
% error coordinate normalized traces (red: LC subsystem, black: anchoring subsystem)
figure;
hold on; 
title("Distance Estimator Trace")
plot((r_polar(202:402,1)-1)./max(r_polar(202:402,1)-1), 'r', 'LineWidth', 2); 
plot(r_polar(202:402,9)./max(r_polar(202:402,9)), 'k', 'LineWidth', 2);
yyaxis right;
plot(([Vielbein(202:402).Metric]./max([Vielbein(202:402).Metric])), 'b', 'LineWidth', 2);
hold off;
