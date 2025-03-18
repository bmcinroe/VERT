function [Vielbein, D, s, r_polar] = computeVert(model, model_params, sim_params, vb_params, window, data)
% Runs the VERT analysis pipeline
% INPUTS:
% model: string indicating either one of the predefined test systems, or set to 'data' if using empirical data
% XX_params: structs containing model, simulation/dataset, and pipeline parameters (see parameter_sweep.m for example)
% window: vector containing initial and final indices of data matrix D for computing VERT. To compute for all data, set: window = [1, length(D)], adjust accordingly to compute over a smaller interval.
% data: a trajectory data matrix D, or empty array if generating data from one of the simulations
% 
% OUTPUTS:
% Vielbein: struct, contains outputs of analysis pipeline
% D: array, data matrix used in VERT computations
% s: array, raw data matrix (same as D unless experimenting with additive noise, etc.)
% r_polar: array, polar coordinate representation (only for multihopf data)


%% Generate data

if strcmp(model, 'cvdp')
    % Run cvdp sim
    cvdpSim();
    r_polar = [];

elseif strcmp(model, 'multihopf')
    % Run multiscale hopf sim'Part
    [s, r_polar, coeff, rlist, tspan] = multiscalehopfSim_fn(model_params, sim_params);

elseif strcmp(model, 'data')
    s = data;
    r_polar = [];
else
    throw('Model not recognized')
end

% Center Data (typically not necessary - leave commented)
% smean = [];
% for i = 1:size(s,2)
%     smean = [smean, mean(s(:,i))];
% end
% 
% D = s - smean;

% Initialize data matrix D, passed to pipeline
D = s;

% Add white Gaussian measurement noise if using
if vb_params.noise == true
    D = awgn(D, vb_params.noise_snr);
end

%% Assign parameters
% unwrap vb_params

% Radius for data partitioning
rho = vb_params.rho;

% Euclidean radius
eps = vb_params.eps;

% Number of nearest neighbors to use for VB analysis
k = vb_params.k;

% Threshold to discard PC's
thresh = vb_params.thresh;

% Time shift for Grassmann/covariance analysis
shift = vb_params.shift;

% Number of trajectories
num_traj = sim_params.numIC * sim_params.numSamp;

% Window size for analysis
tspan = sim_params.tspan;

traj_length = length(tspan);

% Set dt
dt = sim_params.dt;

traj_I = [1];
traj_F = [traj_length];
%dt = 0.01;
for i = 2:num_traj
    traj_I = [traj_I; traj_length*(i-1)+1]; % global indices of trajectory IC's
    traj_F = [traj_F; traj_length*i]; % global indices of trajectory final points
end

%% Partition Samples

% Range of samples to analyze is set by partion_samples

partitions = partition_samples(D, rho, window); % dimension of the input data

%% Find Euclidean Neighborhoods

% Find Euclidean neighborhoods for each sample

[nIndex_euclid, neighborhoods_euclid] = find_neighborhoods(D, eps, partitions, window);

%% kNN Reduction

% Reduce each neighborhood to contain only at most k-Nearest Neighbors

[neighborhoods, nIndex, nSize] = kNN_neighbor_reduction(k, D, nIndex_euclid, neighborhoods_euclid, window);

%% Compute Vielbeins

% Compute the local tangent spaces (VB's) for each sample, the total
% neighborhood variance explained by each VB, and the maximum number of
% non-singular dimensions

% Initialize Vielbein structure for pipeline outputs
Vielbein = struct();

for i = window(1):window(2)
    [nPlane, nExplained, max_tangent_dim, covar] = construct_vielbein(D(i,:), neighborhoods{i}, thresh);

    Vielbein(i).VB = nPlane(:,1:max_tangent_dim);
    Vielbein(i).Explained = nExplained(1:max_tangent_dim);
    Vielbein(i).Dimension = max_tangent_dim;
    % Vielbein(i).Cov = covar; % this is the R^D cov matrix, debug only
    Vielbein(i).Cardinality = size(neighborhoods{i},1);
    Vielbein(i).Neighbors = nIndex{i};

end

VB_dimension = min([Vielbein(window(1):window(2)).Dimension]);

%% Compute VB Grassmann metrics
% Experimental - leave commented

% for i = window(1):window(2) - 1
% 
%     Vielbein(i).PlaneDist = grassDist(Vielbein(i).VB(:,1:VB_dimension), Vielbein(i + 1).VB(:,1:VB_dimension));
% 
% end

% End case
%Vielbein(window(2)).PlaneDist = grassDist(Vielbein(i-1).VB(:,1:VB_dimension), Vielbein(i).VB(:,1:VB_dimension));

%% Compute Anchor Displacement Field

% Compute finite difference estimates of the time derivative for each
% neighborhood, in the observation space

for i = window(1):window(2)

    Vielbein(i).Displacement_i = anchorDisplacements(D, i, traj_I, traj_F, dt); % Displacement for sample i

    Vielbein(i).Displacements = anchorDisplacements(D, nIndex{i}, traj_I, traj_F, dt); % Displacements for neighbors of i

end

%% Project Anchor Displacement Field onto VB's

% Compute finite difference estimates of the time derivative on the VB (tangent space), by
% projecting observation space finite differences onto the VB subspace

for i = window(1):window(2)

    Vielbein(i).VB_Displacements = VBProjection(Vielbein(i).Displacements, Vielbein(i).VB(:,1:Vielbein(i).Dimension));
    Vielbein(i).VB_Displacements = VBProjection(Vielbein(i).Displacements, Vielbein(i).VB(:,1:Vielbein(i).Dimension));
    Vielbein(i).VB_Positions = VBProjection(neighborhoods{i}, Vielbein(i).VB(:,1:Vielbein(i).Dimension));

    Vielbein(i).VB_Position_i = VBProjection(D(i,:), Vielbein(i).VB(:,1:Vielbein(i).Dimension));
    Vielbein(i).VB_Displacement_i = VBProjection(Vielbein(i).Displacement_i, Vielbein(i).VB(:,1:Vielbein(i).Dimension));

    % Compute covariance matrix for projected data
    Vielbein(i).Cov = cov([Vielbein(i).VB_Positions]);


end

%% Compute Template ID Metric

for i = window(1):window(2)

    [metric, beta] = templateMetric(Vielbein(i).VB_Position_i, Vielbein(i).VB_Displacement_i, Vielbein(i).VB_Positions, Vielbein(i).VB_Displacements);

    Vielbein(i).Metric = metric; % Store scalar distance estimate
    Vielbein(i).Jacobian = beta; % Store Jacobian "H(z)", for troubleshooting

end

end

