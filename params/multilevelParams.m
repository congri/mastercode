%params specificially for the multi level model

plt = 1        %plotu_s - u_f regression?
ck = 0         %check?
genTrainingData = 1        %generate new accurate data?
genTestData = 1            %generate new approximate data?

nSurrogates = 3;

nData = 50;                %number of training data
nTest = 50;                %number of test data
nDataFast = 5e3;           %number of fast solver iterations for training data; to construct equally spaced training data
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
trainingDataCov = 1*eye(conductivity.dim);
trainingDataMean = [-2.4948300e+00  -2.3763753e+00  -1.4408964e+00   2.7735040e-01   4.2400721e+00  -4.1972815e+00...
    -6.5863896e+00  -6.0725103e+00  -4.4857574e+00   2.7885775e+00   1.0392218e+00  -4.2395099e-01  -1.7839429e+00...
    -6.4501733e+00   2.7563984e-01  -4.1882800e+00  -6.6131163e+00  -6.0641381e+00  -4.4507977e+00   2.8639372e+00...
    -2.4907763e+00  -2.3980458e+00  -1.2722932e+00   4.7992379e-01   4.2060650e+00];

etaArray = zeros(nSurrogates,bFunOpt.nBasis);