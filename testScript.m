

addpath('./util')
addpath('./MCMC/MCMCsampler')
m = [5 0];
S = eye(size(m, 2));
distribution = @(x) testDistribution(x, m, S);

opts.method = 'MALA'; %randomWalk, nonlocal or langevin
opts.nThermalization = 1e2;
opts.nSamples = 5e4;

%only for random walk
stepWidth = .3;
opts.randomWalk.proposalCov = stepWidth*eye(size(m, 2));

%only for nonlocal
opts.nonlocal.rnd = @(propMean, propCov) mvnrnd(propMean, propCov);
opts.nonlocal.pdf = @(x, propMean, propCov) mvnpdf(x, propMean, propCov);
opts.nonlocal.propMean = zeros(1, size(m,2));
opts.nonlocal.propCov = 10*S;

%only MALA
opts.MALA.stepWidth = 1.2;

startValue = [-10 -2];


[samples, acceptanceRatio] = MCMCsampler(distribution, startValue, opts);
