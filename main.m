%Main script for gradient based MCMC for EM-optimization

clear all;
%load params
addpath('params')
addpath('heatFEM')
addpath(genpath('conductivityField'))
addpath(genpath('MCMC'))
addpath('FEMgradient')
addpath('aux')
addpath('util')
addpath checks
addpath rom
params;
modelType = 'rom';   %which model? reference/multilevel/rom
sv = true;      %save?
dbg = true;     %debug mode? with more output
if(sv)
    disp('   Save data at end of program')
end

%define output function handle giving log U, grad log U and eventually U
if(strcmp(modelType,'reference'))
    
    disp('   Reference solver model')
    nSurrogates = 1;
    outFunction = @(input, conductivity, physical, domain)...
        cont_ref_output(input, conductivity, physical, domain);
    
elseif(strcmp(modelType,'multilevel'))
    
    disp('   Multilevel surrogate model')
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

nIterations = 1;
%preallocate arrays for data recording
prealloc;


%initial value for a in MCMC
if(~exist('aStartInit', 'var'))
    %if there was no thermalization then use regularized mu
    aStartInit = (conductivity.covar_aInv + conductivity.muRegularization*eye(conductivity.dim))...
        *conductivity.covar_a*conductivity.mu_a';
    aStartInit = aStartInit';
end

% If no parallel pool exists
N_Threads = 2;
if isempty(gcp('nocreate'))
    % Create with 2 workers
    parpool('local',N_Threads);
end


%temperature beta
beta = optim.beta(2);

disp('   Started EM optimization...');
t = tic;


%start value for theta
theta = [-6 -6];
thetaArray = theta;
log_sigma_f2 = 0;
K = GPcov(domain.X2, domain.Y2, 'ardSE', [theta, log_sigma_f2]);
Kinv = inv(K);


for nS = 1:nSurrogates
    
    %2 aStart for cold and warm markov chain
    aStart = [aStartInit; aStartInit];
    fixedcovar_a = conductivity.covar_a;
    explorationSteps = 0;
    explorationFactor = 1.1;
    
    
    if(strcmp(modelType,'multilevel'))
        eta = etaArray(nS,:)
        outFunction = @(input, conductivity, physical, domain)...
            surrogateOutput(input, conductivity, physical, domain,...
            etaArray(nS,:), sigmaNoise, bFunH);
    end
    
    
    converged = 0;
    nIterations = 1;
    while(~converged)
        
        %distributions to sample from in M-step
        log_q{1} = @(input) qGP(input, outFunction, conductivity, physical, domain, beta, Kinv);
        log_q{2} = @(input) qGP(input, outFunction, conductivity, physical, domain, beta, Kinv);
        
        
        %Reduce covar_a gradually to make EM more efficient
        if(nIterations <= explorationSteps)
            conductivity.covar_a = explorationFactor^(explorationSteps - nIterations)*fixedcovar_a;
        end
        
        nIterations = nIterations + 1;
        %Total number of MCMC samples in the subsequent step
        if(nIterations < length(opts(1).nSamples))
            nSamplesIteration = opts(1).nSamples(nIterations - 1);
        else
            nSamplesIteration = opts(1).nSamples(end);
        end
        
        %prealloc
        samples = zeros(optim.nSwaps*nSamplesIteration,conductivity.dim);
        
        %M-step: generate samples
        for i = 1:optim.nSwaps
            parfor pp = 1:2
                
                out(pp) = MCMCsampler(log_q{pp}, aStart(pp,:), opts(pp), nSamplesIteration);
                %Refine step width according to acceptance ratio, tune it to about .7
                if(out(pp).acceptance) %has to be nonzero
                    %optim(pp).MCMC.stepWidth = (1/.7)*out(pp).acceptance*optim(pp).MCMC.stepWidth;
                    opts(pp).MALA.stepWidth = (1/.6)*out(pp).acceptance*opts(pp).MALA.stepWidth;
                end
                while(out(pp).acceptance < .03)
                    warning('Acceptance ratio dropped below .03, resample with .3*stepwidth')
                    opts(pp).MALA.stepWidth = .3*opts(pp).MALA.stepWidth;
                    out(pp) = MCMCsampler(log_q{pp}, aStart(pp,:), opts(pp), nSamplesIteration);
                end
                
            end
            
            samples(((i - 1)*nSamplesIteration + 1):(i*nSamplesIteration),:) = out(1).samples;
            if(dbg)
                fprintf(1, '\n       <log(q)>       <log(r)>\n');
                disp([mean(out(1).log_p), mean(out(2).log_p)]);
                fprintf(1, '\n       acc_q       acc_r\n');
                disp([out(1).acceptance, out(2).acceptance]);
                fprintf(1, '\n   step width q    step width r\n');
                disp([opts(1).MALA.stepWidth, opts(2).MALA.stepWidth]);
                fprintf(1, '\n   temperature beta\n');
                disp(beta)
            end
            
            swapLogMetropolis = (1/beta)*out(2).log_pEnd + beta*out(1).log_pEnd - out(1).log_pEnd - out(2).log_pEnd;
            Metropolis = exp(swapLogMetropolis);
            r = rand;
            if(r < Metropolis)
                
                if(dbg)
                    fprintf(1,'\n   State swapping accepted, acceptance ratio:\n');
                    disp(Metropolis);
                end;
                aStart(1,:) = out(2).samples(end,:);
                aStart(2,:) = out(1).samples(end,:);
                optim.betaTrans(2) = optim.betaTrans(2) - 3e-1;
                optim.beta(2) = optim.betaMax(2)*normcdf(optim.betaTrans(2));
            else
                
                if(dbg)
                    fprintf(1,'\n   State swapping rejected, acceptance ratio:\n');
                    disp(Metropolis)
                end
                aStart(1,:) = out(1).samples(end,:);
                aStart(2,:) = out(2).samples(end,:);
                optim.betaTrans(2) = optim.betaTrans(2) + 3e-1;
                optim.beta(2) = optim.betaMax(2)*normcdf(optim.betaTrans(2));
            end
            beta = optim.beta(2);
            log_q{2} = @(input) qGP(input, outFunction, conductivity, physical, domain, beta, Kinv);
            
        end     %samples generated

        



        K_dK = @(lsq) GPcov(domain.X2, domain.Y2, 'ardSE', [lsq, log_sigma_f2]);

        thetaOld = theta;
        decSamples = samples(1:opts(1).nGap:end, :);    %decorrelation
        [theta] = Mstep(decSamples, K_dK, theta)
        thetaArray = [thetaArray; theta]
        K = GPcov(domain.X2, domain.Y2, 'ardSE', [theta, log_sigma_f2]);
        Kinv = inv(K);
        
        fprintf(1,'\n   Number of iterations \n');
        disp(nIterations)
        
        
%         if(norm(theta - thetaOld) < 1e-5)
        if(nIterations > optim.maxIterations)
            converged = 1;
        end
        clear n_diff_Sigma_EM n_diff_mu_EM mu_convergence;
    end
    clear out;
    
%     muOptSurr(nS,:) = conductivity.mu_a;
    
end
clear nS converged swapLogMetropolis fixedcovar_a r Metropolis...
    ans logUMean q_mean grad_q_mean gradLogUMean mu_a_old covar_a_old;





if(strcmp(modelType,'multilevel') && dbg)
    %check true output at last surrogate optimum
    [logUcheck] = cont_ref_output(conductivity.mu_a, conductivity, physical, domain);
    fprintf(1,'\n   log U at surrogate optimum: \n');
    disp(logUcheck);
    clear lambda;
end

computationTime = toc(t); clear t;

if(sv)
    disp('saving...')
    filename = datestr(now,30);
    dir = './data/GP/20';
    mkdir(dir);
    filename = strcat(dir,'/optim',modelType,filename);
    save(filename,'-v7.3');
end

disp('***END***')


