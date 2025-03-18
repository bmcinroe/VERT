clear;
grey = [224 224 224]./255;

%% Initialize model and simulation parameters
% For empirical data, initialize then call 'computeVert' function

% Assign MODEL parameters
model_params = [];
vb_params = [];
sim_params = [];

% dimension of subsystems
model_params.numCycleSys = 1;
model_params.numDecaySys = 8;
model_params.numSys = 2 * model_params.numCycleSys + model_params.numDecaySys;
model_params.dim = 4 * model_params.numCycleSys + 2 * model_params.numDecaySys;

% nonlinearity for polynomial systems
model_params.n = 1;

% frequencies for all subsystems
model_params.param = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

% sweep paramaters
model_params.gamma_rho = 1;
model_params.gamma_theta = 1;
model_params.gamma_d = 1;

% Assign SIM params
sim_params.numIC = 10;
sim_params.numSamp = 5;
sim_params.dt = 0.05;
sim_params.tspan = 0:sim_params.dt:10;
tspan = sim_params.tspan;
sim_params.rho_mean = [2 4]; % lower and upper bounds for radial IC sample
sim_params.numTraj = 10;

% Assign VB params
% Model
model = 'multihopf';
% Radius for data partitioning
vb_params.rho = 1;
% Euclidean radius
vb_params.eps = 0.5;
% Number of nearest neighbors to use for VB analysis
vb_params.k = 30;
% Threshold to discard PC's
vb_params.thresh = 1e-1;
% Time shift for Grassmann/covariance analysis
vb_params.shift = 1;
% Window to compute vert
window = [1, sim_params.numTraj*length(tspan)];
% Noise
vb_params.noise = false;

%% Parameter ranges for \rho and \delta gains
gridSize = 20;
i_range = linspace(1, 3.0, gridSize);
j_range = linspace(1, 3.0, gridSize);
 
[IR, JR] = meshgrid(i_range, j_range);

%% Generate baseline submanifold data
% Used to assess accuracy of VERT outputs for Hopf model system

% Generate cycle ground truth dataset

attractor_points = 1000;

rho_1 = ones(1,attractor_points).';
rho_2 = rho_1;

theta_1 = linspace(0,2*pi,attractor_points).';
theta_2 = theta_1;

s_polar_groundtruth = [rho_1, theta_1, rho_2, theta_2];

for i = 1:model_params.numDecaySys
    s_polar_groundtruth = [s_polar_groundtruth, zeros(1, attractor_points).', zeros(1, attractor_points).'];
end

s_groundtruth = zeros(size(s_polar_groundtruth));

% Iterate through each plane and tranform to cartesian coordinates
for i = 1:2:(size(s_polar_groundtruth, 2))
    [x, y] = pol2cart(s_polar_groundtruth(:,i+1), s_polar_groundtruth(:,i));
    s_groundtruth(:,i) = x;
    s_groundtruth(:,i+1) = y;
end

% Generate torus ground truth dataset
attractor_points_torus = 30;

theta_torus = linspace(0,2*pi,attractor_points_torus).';

[theta1_torus, theta2_torus] = meshgrid(theta_torus, theta_torus);

torus_polar = [];
for i = 1:attractor_points_torus
    for j = 1:attractor_points_torus
        torus_polar = [torus_polar; 1, theta1_torus(i,j), 1, theta2_torus(i,j)];
    end
end

[x1, y1] = pol2cart(torus_polar(:,2), torus_polar(:,1));
[x2, y2] = pol2cart(torus_polar(:,4), torus_polar(:,3));

torus_groundtruth = zeros(length(x1), model_params.dim);
for i = 1:size(torus_groundtruth, 1)
    torus_groundtruth(i,:) = [x1(i) y1(i) x2(i) y2(i) zeros(1, model_params.dim - 4)];
end

baselines.cycle = s_groundtruth;
baselines.torus = torus_groundtruth;

%% Parameter Sweep Loop

% Number of trials for each parameter combination
numTrials = 1;

% Cache for parameter combos that produce error (for troubleshooting, should return empty)
failed_trials = {};

% Store metric and trajectory data for all trials
vb_metric = cell(length(i_range),length(j_range));
vb_trajectory = cell(length(i_range),length(j_range));

tic
% Replace 'parfor'-> 'for' to convert to standard for-loop
% Non-iterable parameters can be assigned outside if using sequential loop
parfor ii = 1:numel(IR)

    trial = 1;

    vb_metric_trial = [];

    test_trajectories = {numTrials};

    while trial <= numTrials
        % Assign parameters for each worker if using parallel loop
        % Assign MODEL parameters
        model_params = [];
        vb_params = [];
        sim_params = [];

        % dimension of subsystems
        model_params.numCycleSys = 1;
        model_params.numDecaySys = 8;
        model_params.numSys = 2 * model_params.numCycleSys + model_params.numDecaySys;
        model_params.dim = 4 * model_params.numCycleSys + 2 * model_params.numDecaySys;

        % nonlinearity for polynomial systems
        model_params.n = 1;

        % frequencies for all subsystems
        model_params.param = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

        % sweep paramaters
        model_params.gamma_rho = JR(ii);
        model_params.gamma_theta = 1;
        model_params.gamma_d = IR(ii);

        % Assign SIM params
        sim_params.numIC = 10;
        sim_params.numSamp = 5;
        sim_params.dt = 0.05;
        sim_params.tspan = 0:sim_params.dt:10;
        tspan = sim_params.tspan;
        sim_params.rho_mean = [2 4]; % lower and upper bounds for radial IC sample
        sim_params.numTraj = 10; % number of trajectories to analyze

        % Assign VB params
        % Model
        model = 'multihopf';
        % Radius for data partitioning
        vb_params.rho = 1;
        % Euclidean radius
        vb_params.eps = 0.5;
        % Number of nearest neighbors to use for VB analysis
        vb_params.k = 30;
        % Threshold to discard PC's
        vb_params.thresh = 1e-1;
        % Time shift for Grassmann/covariance analysis
        vb_params.shift = 1;
        % Window to compute vert
        window = [1, sim_params.numTraj*length(tspan)];
        % Noise
        vb_params.noise = false;
        vb_params.noise_snr = 60;

        try
            [Vielbein, D, s, r_polar] = computeVert(model, model_params, sim_params, vb_params, window, []);
        catch
            fprintf(strcat('Analysis failed for trial ', string(ii)));
            failed_trials = [failed_trials; [IR(ii),JR(ii)]];
        end

        % Store VB metric vector
        vb_metric_trial = [vb_metric_trial; [Vielbein.Metric]];

        % Store tested trajectory
        test_trajectories{trial,1} = D(window(1):window(2), :);
        
        % Increment trial counter
        trial = trial + 1;
    end

    % Add raw trajectory data and VERT output
    vb_metric{ii} = vb_metric_trial;
    vb_trajectory{ii} = test_trajectories;

end
toc