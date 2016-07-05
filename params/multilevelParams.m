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
trainingDataMean = [-1.8271070e+00  -3.3049167e+00  -1.0888289e+00   6.9051042e-01   2.9253278e+00  -4.7117301e+00  -6.9981758e+00  -5.7406818e+00  -1.0495619e+00   1.5646482e+00  -2.1002816e+00  -3.6223970e+00  -5.9550131e+00  -2.3994230e+00  -1.4386804e-01  -4.7146378e+00  -6.8978550e+00  -5.7035896e+00  -1.1187027e+00   1.4617227e+00  -1.8094909e+00  -3.4121527e+00  -1.0525894e+00   9.5777496e-01   2.8658767e+00];

etaArray = zeros(nSurrogates,bFunOpt.nBasis);