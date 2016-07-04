

addpath('./util')
addpath('./MCMC/MCMCsampler')
m = [5 0];
S = eye(size(m, 2));
S(2,2) = 9;
distribution = @(x) testDistribution(x, m, S);
% gradp = @(x) testDistributionGrad(x, m, S);

opts.method = 'MALA'; %randomWalk, nonlocal or langevin
opts.nThermalization = 0e2;
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

startValue = [-100 -20];


[out] = MCMCsampler(distribution, startValue, opts);

%HMC
% options(1) = 0;             %display
% options(5) = 1;             %momentum persistence
% options(7) = 2;             %leapfrog steps
% options(9) = 0;             %gradient check
% options(14) = 5e4;          %nSamples
% options(15) = 0;            %nTherm
% options(17) = 1;           %persistence?
% options(18) = .2/options(7); %step size in leapfrogs
% [samplesHMC, energies, diagn] = hmc(distribution, startValue, options, gradp);
