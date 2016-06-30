%train a multilevel model as a surrogate

multilevelParams;

%Random seed
rng('shuffle')



%accurate data
if(genTrainingData)
    
    %generate a dataset of a's
    inputData = mvnrnd(trainingDataMean,trainingDataCov,nDataFast);
    
    %save input data
    mkdir('./data/multilevel')
    namestr = strcat('./data/multilevel/','input');
    save(namestr, 'inputData','-v7.3');
    
    %specify mesh size before execution of params script!
    clear domain;
    domain.nElements = nFast; %number of finite elements
    paramsNewDomain;
    [~, approx_u] =...
        generateData(inputData, physical, domain, conductivity, nDataFast, 'approxData');
    
    [approx_u, inputData] = pickData(approx_u, inputData, nData, upperFraction);
    
    clear domain;
    domain.nElements = nSlow; %number of finite elements
    paramsNewDomain;
    [accurateExp, accurate_u] =...
        generateData(inputData, physical, domain, conductivity, nData, 'accurateData');
    
else
    temp = load('./data/multilevel/accurateData.mat');
    accurateExp = temp.logU;
    accurate_u = temp.transformedU;
    temp = load('./data/multilevel/approxData.mat');
%     approxExp = temp.logU;
    approx_u = temp.transformedU;
    temp = load('./data/multilevel/input.mat');
    inputData = temp.inputData;
    [approx_u, inputData] = pickData(approx_u, inputData, nData, upperFraction);

    clear temp;
    
    %set to fast solver mesh size
    clear domain;
    domain.nElements = nFast; %number of finite elements
    paramsNewDomain;
end

%sort data
[approx_u, index] = sort(approx_u);
accurate_u = accurate_u(index);
% [approxExp, index] = sort(approxExp);
% accurateExp = accurateExp(index);
clear index;

if(strcmp(bFunOpt.type, 'rbf') || strcmp(bFunOpt.type, 'rbf+trend'))
    meanApprox_u = mean(approx_u);
    widthApprox_u = max(approx_u) - min(approx_u);
    bFunOpt.rbfLocations = linspace(meanApprox_u - .2*widthApprox_u, meanApprox_u + .2*widthApprox_u, bFunOpt.nBasis); %only for rbf models
    bFunOpt.rbfStd = (max(approx_u) - min(approx_u))/bFunOpt.nBasis;
end

%basis function handle
bFunH = @(inp)basisFunctions(inp, bFunOpt);
%design matrix
Phi = bFunH(approx_u)';

sigmaNoise = optimNoise(accurate_u, Phi, priorMean, priorCov, mu_sigma)
% sigmaNoise = 100

%posterior params on basis function coeffs
[postMean, postCovInv] = postParams(priorMean, priorCovInv, sigmaNoise, Phi, accurate_u);
postCov = inv(postCovInv);
%model properly set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ck)
    eta = mvnrnd(postMean, postCov)
    %perform check: surrogate output vs. true output
    nCheck = 50;
    aSamples = mvnrnd(zeros(1,conductivity.dim),9*eye(conductivity.dim), nCheck);
    for i = 1:nCheck
        i
        [lambda] = computeLambda(aSamples(i,:), domain, conductivity.lambdaCutoff);
        [surrLogUi, surrGradLogUi] = surrogateOutput(lambda, aSamples(i,:), conductivity, physical, domain,...
     eta, sigmaNoise, bFunH);
    surrLogU(i) = surrLogUi;
    end
    clear domain;
    domain.nElements = nSlow; %number of finite elements
    paramsNewDomain;
    for i = 1:nCheck
        i
        [lambda] = computeLambda(aSamples(i,:), domain, conductivity.lambdaCutoff);
        [logU, gradLogU] = cont_ref_output(lambda, aSamples(i,:), conductivity, physical, domain);
        trueExp(i) = logU;
    end
    %set domain back to coarse
    clear domain;
    domain.nElements = nFast;
    paramsNewDomain;

    [surrLogU' trueExp']
end
clear ck;


if(genTestData)
    inputTest = mvnrnd(trainingDataMean,trainingDataCov,nTest);
    %specify mesh size before execution of params script!
    clear domain;
    domain.nElements = nSlow; %number of finite elements
    paramsNewDomain;
    [accurateExpTest, accurate_uTest] =...
        generateData(inputTest, physical, domain, conductivity, nTest, 'accurateDataTest');
    clear domain;
    domain.nElements = nFast; %number of finite elements
    paramsNewDomain;
    [approxExpTest, approx_uTest] =...
        generateData(inputTest, physical, domain, conductivity, nTest, 'approxDataTest');
else
    temp = load('./data/multilevel/accurateDataTest.mat');
    accurateExpTest = temp.logU;
    accurate_uTest = temp.transformedU;
    temp = load('./data/multilevel/approxDataTest.mat');
    approxExpTest = temp.logU;
    approx_uTest = temp.transformedU;
    clear temp;
end
%sort
[approx_uTest, index] = sort(approx_uTest);
accurate_uTest = accurate_uTest(index);
[approxExpTest, index] = sort(approxExpTest);
accurateExpTest = accurateExpTest(index);
clear index;

approxInf = linspace(min(approx_u) - .3*abs(min(approx_u)), max(approx_u) + .3*abs(min(approx_u)), 500);
% [predMean, predVar] = predParams(postMean, postCov, approxInf', sigmaNoise, bFunH);



if(plt)
    
    figure
%     erp = shadedErrorBar(approxInf,predMean,predVar)
    hold on
    pac = plot(approx_uTest,accurate_uTest,'x');
    pap = plot(approx_u,accurate_u,'x');
    for i = 1:4
        eta = mvnrnd(postMean, postCov);
        samp = eta*bFunH(approxInf);
        plot(approxInf, samp,'k')
    end
    hold off
    box on
    clear eta approxInf samp;
    
end
clear plt;



