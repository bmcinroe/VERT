%% Initialize IC's

IdenticalPopulationMeans = true;

numIC = 20; % Number of IC means
numSamp = 1; % Number of samples 

IC_list = cell(numIC, 1);

var = 0.2; % Variance of distribution for IC's

% Populate Initial Condition Set
for i = 1:numIC
    
    IC_set = zeros(numSamp,4);
    state_mean = [normrnd(0.5, var), normrnd(0.3, var), normrnd(0.5, var), normrnd(0.1, var)]';
    
    for j = 1:numSamp
        
        IC = state_mean + [normrnd(0.0, var), normrnd(0.0, var), normrnd(0.0, var), normrnd(0.0, var)]';
        IC_set(j,:) = IC;
    end
    IC_list{i} = IC_set;
end


%% Integrate Dynamics

tspan = 0:0.1:20;

s = [];

for i = 1:numIC

    for j = 1:numSamp

        r0 = IC_list{i}(j,:); 
        [t,r_i] = ode45('cvdp',tspan,r0,[]);
        s = [s; r_i];

    end

end

[coeff,score,latent,tsquared,explained,mu] = pca(s);

%plot3(s * coeff(:,1), s * coeff(:,2), s * coeff(:,3), '.');