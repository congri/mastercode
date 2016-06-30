%params specificially for the multi level model

plt = 1        %plotu_s - u_f regression?
ck = 0         %check?
genTrainingData = 1        %generate new accurate data?
genTestData = 1            %generate new approximate data?

nSurrogates = 3;

nData = 50;                %number of training data
nTest = 50;                %number of test data
nDataFast = 2e3;           %number of fast solver iterations for training data; to construct equally spaced training data
upperFraction = .5;       %only take upper fraction of fast solver data to build surrogate, as we mainly want to
                           %predict high values of u_s

bFunOpt.type = 'rbf+trend';
bFunOpt.nBasis = 4;
%mean of exponential prior on sigma, must be > 0
mu_sigma = 1e100;
%linear model; posterior mean and cov of model params
priorMean = zeros(bFunOpt.nBasis,1); %assuming that accurate = approximate
priorCov = 1e6*eye(bFunOpt.nBasis); %choose a broad prior
priorCovInv = inv(priorCov);

%nElements contains number of elements of top level code in 1 and fastest
%code in last element
nFast = 225;
nSlow = domain.nElements;

%training data params
trainingDataCov = 36*eye(conductivity.dim);
trainingDataMean = zeros(1,conductivity.dim);

etaArray = zeros(nSurrogates,bFunOpt.nBasis);