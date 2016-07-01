
m = [1 0];
S = eye(2);
distribution = @(x) testDistribution(x, m, S);

opts.method = 'nonlocal'; %randomWalk, nonlocal or langevin
opts.nThermalization = 1e2;
opts.nSamples = 5e4;

%only for random walk

%only for nonlocal
opts.nonlocalrnd = @() mvnrnd(zeros(1, size(m, 2)), S);
opts.nonlocalpdf = @(x) mvnpdf(x, zeros(1, size(m, 2)), S);

startValue = 4*rand(1, 2) - 2;


[samples, acceptanceRatio] = MCMCsampler(distribution, startValue, opts);
