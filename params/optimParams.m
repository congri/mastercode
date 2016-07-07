%params for expectation maximization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optim.nSwaps = 5; %how often propose to swap markov chains

%maximum iterations in opt step
optim.maxIterations = 340;
%Convergence criteria
optim.muDifferenceThreshold = 1e-20;
optim.covarDifferenceThreshold = 1e-4;
optim.betaTrans = Inf;
optim.beta = 1;

%second, tempered Markov chain
optim.betaTrans = [Inf 2];
optim.betaMax = [1 1];
optim.beta = [1 optim.betaMax*normcdf(optim.betaTrans(2))];

%MCMC settings
opts.method = 'MALA';              %proposal type: randomWalk, nonlocal or MALA
opts.nThermalization = 0;      %thermalization steps
opts.nSamples = 200;             %number of samples

% only for MALA
opts.MALA.stepWidth = 3e-1;       %step size parameter
opts = repmat(opts,2,1);
opts(2).MALA.stepWidth = opts(1).MALA.stepWidth;

