%params specificially for the multi level model

plt = 0        %plotu_s - u_f regression?
ck = 0         %check?
genTrainingData = 0        %generate new accurate data?
genTestData = 0            %generate new approximate data?

nSurrogates = 3;

nData = 50;                %number of training data
nTest = 50;                %number of test data
nDataFast = 2e3;           %number of fast solver iterations for training data; to construct equally spaced training data
upperFraction = 1;       %only take upper fraction of fast solver data to build surrogate, as we mainly want to
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
trainingDataCov = 1*eye(conductivity.dim);
trainingDataMean = [-6.1761148e-01  -4.2258421e+00  -3.6725914e+00   1.3968711e-01   4.9445302e+00  -3.8505849e+00...
    -7.1101403e+00  -6.7356088e+00  -5.0692089e+00   2.6493758e+00  -2.0053591e+00  -4.1576103e+00  -4.4035804e+00...
    -6.6025675e+00  -2.7257838e-01  -3.8531875e+00  -7.1306576e+00  -6.9422608e+00  -5.1534670e+00   2.9569076e+00...
    -5.9642452e-01  -4.2334568e+00  -3.7238121e+00   4.5064811e-02   4.8454648e+00];

etaArray = zeros(nSurrogates,bFunOpt.nBasis);