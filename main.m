%Main script for gradient based MCMC for EM-optimization


clear all;
%load params
addpath('params')
addpath('heatFEM')
addpath('conductivityField')
addpath('MCMC')
addpath('FEMgradient')
addpath('aux')
addpath('util')
params;
modelType = 'multilevel';   %which model? reference/multilevel
sv = false      %save?

%define output function handle giving log U, grad log U and eventually U
if(strcmp(modelType,'reference'))
    
    disp('Reference solver model')
    nSurrogates = 1;
    outFunction = @(lambda, input, conductivity, physical, domain)...
        cont_ref_output(lambda, input, conductivity, physical, domain);
    
elseif(strcmp(modelType,'multilevel'))
    
    disp('Multilevel surrogate model')
    multilevelParams;
    addpath('multilevel')
    addpath('multilevel/linearModel')
    %train the model
    trainMultilevel;
    etaArray = mvnrnd(postMean, postCov, nSurrogates);
    
    error('model trained') %to stop after model training

    
else
    
    error('Unknown model type')
    
end




disp('Started EM optimization...');
t = tic;

nIterations = 1;
%preallocate arrays for data recording
prealloc;


%initial value for a in MCMC
if(~exist('aStart', 'var'))
    %if there was no thermalization then use regularized mu
    aStart = (conductivity.covar_aInv + conductivity.muRegularization*eye(conductivity.dim))*conductivity.covar_a*conductivity.mu_a';
    aStart = aStart';
end
%2 aStart for cold and warm markov chain
aStart = [aStart; aStart];
fixedcovar_a = conductivity.covar_a;
explorationSteps = 0;
explorationFactor = 1.1;

% If no parallel pool exists
N_Threads = 2;
if isempty(gcp('nocreate'))
    % Create with 2 workers
    parpool('local',N_Threads);
end

for nS = 1:nSurrogates

    if(strcmp(modelType,'multilevel'))
        etaArray(nS,:)
        outFunction = @(lambda, input, conductivity, physical, domain)...
            surrogateOutput(lambda, input, conductivity, physical, domain,...
         etaArray(nS,:), sigmaNoise, bFunH);
    end
    
    
    converged = 0;
    nIterations = 1;
    while(~converged)
        %Reduce covar_a gradually to make EM more efficient
        if(nIterations <= explorationSteps)
            conductivity.covar_a = explorationFactor^(explorationSteps - nIterations)*fixedcovar_a;
        end
        clear fixedcovar_a;

        nIterations = nIterations + 1;
        %Total number of MCMC samples in the subsequent step
        if(nIterations < length(optim(1).MCMC.nSamples))
            nSamplesIteration = optim(1).MCMC.nSamples(nIterations - 1);
        else
            nSamplesIteration = optim(1).MCMC.nSamples(end);
        end

        %prealloc
        samples = zeros(optim(1).nSwaps*nSamplesIteration,conductivity.dim);
        logUMean = 0;
        q_mean = 0;
        grad_q_mean = zeros(1,conductivity.dim);
        gradLogUMean = zeros(1,conductivity.dim);

        %M-step
        for i = 1:optim(1).nSwaps
            parfor pp = 1:2
                out(pp) = a_samples(outFunction, physical, conductivity, optim(pp),...
                    domain, nSamplesIteration, aStart(pp,:));
                optim(pp) = out(pp).optim;
                %Refine step width according to acceptance ratio, tune it to about .7
                if(out(pp).acceptance) %has to be nonzero
                    optim(pp).MCMC.stepWidth = (1/.7)*out(pp).acceptance*optim(pp).MCMC.stepWidth;
                end
                while(out(pp).acceptance < .03)
                warning('Acceptance ratio dropped below .03, resample with .3*stepwidth')
                optim(pp).MCMC.stepWidth = .3*optim(pp).MCMC.stepWidth;
                out(pp) = a_samples(outFunction, physical, conductivity, optim(pp),...
                    domain, nSamplesIteration, aStart(pp,:));
                end
                chain = pp
                chainAcceptance = out(pp).acceptance
            end
            
            samples(((i - 1)*nSamplesIteration + 1):(i*nSamplesIteration),:) = out(1).samples;
            logUMean = ((i - 1)/i)*logUMean + (1/i)*out(1).statistics.logUmean
            q_mean = ((i - 1)/i)*q_mean + (1/i)*out(1).statistics.q_mean;
            grad_q_mean = ((i - 1)/i)*grad_q_mean + (1/i)*out(1).statistics.grad_q_mean;
            gradLogUMean = ((i - 1)/i)*gradLogUMean + (1/i)*out(1).statistics.gradLogUMean;
            stepWidth = optim(1).MCMC.stepWidth
            stepWidthTempered = optim(2).MCMC.stepWidth

            

            swapLogMetropolis = (optim(1).beta - optim(2).beta)*(out(2).logU - out(1).logU + ...
                out(2).pExponent - out(1).pExponent)
            Metropolis = exp(swapLogMetropolis); clear swapLogMetropolis;
            r = rand;
            if(r < Metropolis)
               disp('State swapping accepted') 
               aStart(1,:) = out(2).samples(end,:);
               aStart(2,:) = out(1).samples(end,:);
               optim(2).betaTrans = optim(2).betaTrans - 1e-1;
               optim(2).beta = optim(2).betaMax*normcdf(optim(2).betaTrans);
            else
               disp('State swapping rejected')
               aStart(1,:) = out(1).samples(end,:);
               aStart(2,:) = out(2).samples(end,:);
               optim(2).betaTrans = optim(2).betaTrans + 1e-1;
               optim(2).beta = optim(2).betaMax*normcdf(optim(2).betaTrans);
            end
            beta = optim(2).beta
            clear r Metropolis;
            
        end
        clear i;

        mu_a_old = conductivity.mu_a;

        %Regularization
        conductivity.mu_a = (conductivity.covar_aInv + ...
            conductivity.muRegularization*eye(conductivity.dim))\(conductivity.covar_aInv*(mean(samples,1))');
        conductivity.mu_a = conductivity.mu_a';
        %only keep a single component of mu this time
%         mu = conductivity.mu_a(component)
%         conductivity.mu_a = mu_a_old;
%         conductivity.mu_a(component) = mu; clear mu;
        conductivity.mu_a
        clear ans;

        covar_a_old = conductivity.covar_a;

        muArray(nIterations,:) = conductivity.mu_a;
        logUmeanArray(nIterations) = logUMean; clear logUMean;
        q_meanArray(nIterations) = q_mean; clear q_mean;
        grad_q_meanArray(nIterations,:) = grad_q_mean; clear grad_q_mean;
        gradLogUMeanArray(nIterations,:) = gradLogUMean; clear gradLogUMean;

        %convergence condition
        n_diff_mu_EM = norm(mu_a_old - conductivity.mu_a); clear mu_a_old;
        mu_convergence = n_diff_mu_EM/conductivity.dim
        nIterations
        n_diff_Sigma_EM = norm(covar_a_old - conductivity.covar_a);
        clear covar_a_old;
        muDiffArray(nIterations) = n_diff_mu_EM;
        if(((mu_convergence < optim(1).muDifferenceThreshold &&...
                n_diff_Sigma_EM/(conductivity.dim^2) < optim(1).covarDifferenceThreshold) && nIterations > explorationSteps)...
                || nIterations > optim(1).maxIterations)
            converged = 1;
        end
        clear n_diff_Sigma_EM n_diff_mu_EM mu_convergence;
    end
    clear out;

    muOptSurr(nS,:) = conductivity.mu_a;
    
end
clear nS converged;

computationTime = toc(t); clear t;

if(strcmp(modelType,'multilevel'))
    %check true output at last surrogate optimum
    lambda = computeLambda(conductivity.mu_a, domain, conductivity.lambdaCutoff);
    [logUcheck] = cont_ref_output(lambda, conductivity.mu_a, conductivity, physical, domain)
    clear lambda;
end

%clear auxiliary variables for clarity
if(sv)
    disp('saving...')
    filename = datestr(now,30);
    dir = './data/dim9/comp6';
    mkdir(dir);
    filename = strcat(dir,'/optim',modelType,filename);
    save(filename,'-v7.3');
end

disp('***END***')


