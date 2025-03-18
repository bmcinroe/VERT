function rhs = multiscalehopf_fn(t,r,dummy,model_params)
% ODE for multiscale coupled Hopf Oscillator system 

% dimension of subsystems
numCycleSys = model_params.numCycleSys;
numDecaySys = model_params.numDecaySys;
numSys = model_params.numSys;
dim = model_params.dim;

% subsystem gains
gamma_rho = model_params.gamma_rho;
gamma_theta = model_params.gamma_theta;
gamma_d = model_params.gamma_d;

% nonlinearity for polynomial systems
n = model_params.n;

% frequencies for all subsystems
param = model_params.param;

% codimension of template relative to the anchor
codim = dim - (numCycleSys * 4);

% radius of limit cycle
radius = 1.0;

rhs = [];

% Generate rhs for limit cycle subsystems
for i = 1:numCycleSys
    ind1 = 2*i - 1;
    ind2 = 2*(i + 1) - 1;
    rhs = [rhs; gamma_rho * (radius - r(ind1))*r(ind1); param(i) + gamma_theta * sin(r(2*(i+1)) - r(2*i))]; % add first module
    rhs = [rhs; gamma_rho * (radius - r(ind2))*r(ind2); param(i+1) + gamma_theta * sin(r(2*i) - r(2*(i+1)))]; % add second module
end

% Generate rhs for decaying subsystems
for i = (numCycleSys + 2):numSys
    ind = 2*i - 1;
    rhs = [rhs; -gamma_d * r(ind)^n; param(i)];
end

% Uncomment if frequencies are state dependent

% beta = param(2) * sin(r(2))^9;
% gamma = param(4) * sin(r(2))^9;

end