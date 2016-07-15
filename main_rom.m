%Main file for reduced order model

clear all
addpath rom
addpath params
addpath heatFEM
addpath util
addpath(genpath('conductivityField'))
addpath(genpath('MCMC'))




nFast = 25;
nSlow = 1600;
domain.nElements = nSlow;

n_f = [sqrt(nSlow) sqrt(nSlow)];
n_c = [sqrt(nFast) sqrt(nFast)];
[E] = get_coarse_el(n_f, n_c, 1:nSlow);

%load params/ precomputation
params;

%data
nData = 4;
input = zeros(nData, nSlow);
%T_fine is n_natNodes x nData
nnodes_fine = (sqrt(nSlow) + 1)^2;
T_fine = zeros(nnodes_fine, nData);
for i = 1:nData
    input(i, :) = mvnrnd(zeros(1, nSlow), 9*eye(nSlow));
    [T_fine(:, i)] = romOutput(input(i, :), conductivity, physical, domain);
end

%Switch to fast solver
%specify mesh size before execution of params script!
clear domain;
domain.nElements = nFast; %number of finite elements
paramsNewDomain;
nnodes_coarse = (sqrt(nFast) + 1)^2;


%EM start values
%X might not be correct, just testing here
% xk = input(1, :);
% xk = reshape(xk, nSlow/nFast, nFast);

%inputs belonging to same coarse element are in the same column of xk
xk = zeros(nSlow/nFast, nFast);
for i = 1:nFast
    xk(:,i) = input(1, E == i)';
end

theta_c.theta = [1; 1];
theta_c.varInv = 4;
theta_cf.mu = zeros(nnodes_fine, 1);
theta_cf.W = eye(nnodes_fine, nnodes_coarse);
theta_cf.Sinv = eye(nnodes_fine);
%design matrix Phi
Phi = [mean(xk); var(xk)];

%start at mean
X = theta_c.theta'*Phi;



% If no parallel pool exists
N_Threads = 2;
if isempty(gcp('nocreate'))
    % Create with 2 workers
    parpool('local',N_Threads);
end

%separate q's for every data point
q = cell(nData, 1);
for i = 1:nData
    q{i} = @(Xi) LBdist(Xi, T_fine(:, i), Phi, theta_c, theta_cf, physical, domain, conductivity);
end

clear opts;
opts.method = 'randomWalk';
stepWidth = 8e-3;
opts.randomWalk.proposalCov = stepWidth*eye(size(X, 2));
opts.nThermalization = 0;
opts.nSamples = 1000;

parfor pp = 1:nData
    out(pp) = MCMCsampler(q{pp}, X, opts);
end


