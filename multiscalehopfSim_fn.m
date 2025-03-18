function [s, r_polar, coeff, rlist, tspan] = multiscalehopfSim_fn(model_params, sim_params)
% return cartesian state, polar state, and global PCA coefficients


%% Generalized Hopf System
% Initialize data

numCycleSys = model_params.numCycleSys;
numDecaySys = model_params.numDecaySys;
numSys = model_params.numSys;
dim = model_params.dim;

numIC = sim_params.numIC;
numSamp = sim_params.numSamp;

tspan = sim_params.tspan;

rho_mean = sim_params.rho_mean;

rho_mean_lower = rho_mean(1);
rho_mean_upper = rho_mean(2);

% Uniformly sample IC means on the unit disk for each subsystem

IdenticalPopulationMeans = true; % All oscillators in subpopulation (LC / decay) have same means

IC_list = cell(numIC, 1);

for i = 1:numIC
    
    rho_0 = zeros(numSys,1);
    theta_0 = zeros(numSys,1);
    
    % populate IC means for LC subsystems
    rho_mean = unifrnd(rho_mean_lower, rho_mean_upper);
    theta_mean = unifrnd(0, 2*pi);
    for j = 1:numCycleSys
        rho_0(j) = rho_mean;
        theta_0(j) = theta_mean;
    end
    
    % populate IC means for decaying subsystems
    rho_mean = unifrnd(rho_mean_lower, rho_mean_upper);
    theta_mean = unifrnd(0, 2*pi);
    for j = (numCycleSys + 1):numSys
        rho_0(j) = rho_mean;
        theta_0(j) = theta_mean;
    end
    
    IC = [rho_0 theta_0]';
    
    IC_list{i} = IC(:);
    
end

rlist = {};

% Populate cells with initial conditions
for k = 1:length(IC_list)
    for i = 1:numSamp
        ri = zeros(dim,1);
        for j = 1:size(ri,1)
            ri(j) = IC_list{k}(j) + normrnd(0,0.3); % base sig^2 = 0.3
        end
        rlist = [rlist; ri];
    end
end

% Integrate for each initial condition, save trajectory in r
r = [];
for i = 1:length(rlist)
    [ti,ri] = ode45('multiscalehopf_fn',tspan,rlist{i},[],model_params);
    r = [r; ri];
end

r_polar = r;

% Transform trajectories to cartesian coordinates

% initalize cartesian trajectories
s = zeros(size(r));

% Iterate through each plane and tranform to cartesian coordinates
for i = 1:2:(dim)
    [x, y] = pol2cart(r(:,i+1), r(:,i));
    s(:,i) = x;
    s(:,i+1) = y;
end



% PCA projection of global attractor
[coeff,score,latent,tsquared,explained,mu] = pca(s);

%plot3(s * coeff(:,1), s * coeff(:,2), s * coeff(:,3), '.');

end