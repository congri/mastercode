%params for expectation maximization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MCMC step variance; initial value for local random walk
optim.samplingMethod = 'langevin';
optim.MCMC.stepWidth = 3e-3;
%Number of MCMC samples in the first iterations; every following iteration
%has as much samples as the last element of nSamples
% optim.MCMC.nSamples = 5*ceil(1.1.^(1:20));
% optim.MCMC.nTherm = 10; %MCMC thermalization; same in each step
optim.MCMC.nSamples = [100];%[10 20 30 40 50 75 100 200 300 400];
optim.MCMC.nTherm = 0;
optim.MCMC.nGap = 1;  %MCMC gap, same in each step
optim.MCMC.runTimeLimit = 1000000; %Runtime limit for single sampling step
optim.nSwaps = 5; %how often propose to swap markov chains

%maximum iterations in opt step
optim.maxIterations = 70;
%Convergence criteria
optim.muDifferenceThreshold = 1e-20;
optim.covarDifferenceThreshold = 1e-4;
optim.betaTrans = Inf;
optim.beta = 1;

%nChains optim structures for nChains parallel Markov chains
nChains = 2;
optim = repmat(optim,2,1);
%second, tempered Markov chain
optim(2).MCMC.stepWidth = 10*optim(1).MCMC.stepWidth;
optim(2).betaTrans = -2;
optim(2).betaMax = 1;
optim(2).beta = optim(2).betaMax*normcdf(optim(2).betaTrans);

%MCMC settings
opts.method = 'MALA';              %proposal type: randomWalk, nonlocal or MALA
opts.nThermalization = 10;      %thermalization steps
opts.nSamples = 100;             %number of samples

% only for MALA
opts.MALA.stepWidth = 3e-2;       %step size parameter
opts = repmat(opts,2,1);
opts(2).MALA.stepWidth = 10*opts(1).MALA.stepWidth;

