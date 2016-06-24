function [out] = a_samples(outFunction, physical, conductivity, optim, domain, nSamplesIteration, a)
%Compute samples of q^beta(a) = (U(a)*p(a|theta))^beta in EM algorithm,
%where beta might be an auxiliary temperature for multi-modality mixing

rng('shuffle');
%measure runtime
tic;

sampleSize = ceil(nSamplesIteration/optim.MCMC.nGap);
samples = NaN*zeros(sampleSize,numel(conductivity.mu_a));

%zero mean and identity covariance for random z in langevin proposal
zMean = zeros(1, numel(conductivity.mu_a));
zCovariance = eye(numel(conductivity.mu_a));
invProposalCov = inv(optim.MCMC.stepWidth^2*zCovariance);

%Values of heat conductivities on finite elements
lambda = computeLambda(a, domain, conductivity.lambdaCutoff);

[logU, gradLogU] = outFunction(lambda, a, conductivity, physical, domain);
aDiff = (a - conductivity.mu_a);
pExponent = -.5*optim.beta*aDiff*conductivity.covar_aInv*aDiff';
%assert(q ~= 0, 'Error: MCMC weight q dropped to 0 numerically')
%assert(outFunctionOut.output ~= 0, 'Error: MCMC weight q dropped to 0 numerically')
assert(~isinf(logU), 'Error: infinite output exponent: zero probability')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Thermalization
accepted_therm = 0;
accThermRatio = 0;
statistics.logUmean = 0;
statistics.q_mean = 0;
grad_q = -(conductivity.covar_a\((a - conductivity.mu_a)'))' + gradLogU;
statistics.grad_q_mean = grad_q;
statistics.gradLogUMean = gradLogU;
for i = 1:optim.MCMC.nTherm
    if(strcmp(optim.samplingMethod,'langevin'))
        %Langevin proposal 15/02/16
        proposalMean = a + .5*optim.beta*optim.MCMC.stepWidth^2*grad_q;
        aTemp = proposalMean + optim.MCMC.stepWidth*mvnrnd(zMean, zCovariance);

        lambda = computeLambda(aTemp, domain, conductivity.lambdaCutoff);

        [logUprop, gradLogUprop] = outFunction(lambda, aTemp, conductivity, physical, domain);
        grad_qprop = -(conductivity.covar_a\((aTemp - conductivity.mu_a)'))' + gradLogUprop;

        proposalExponent = -.5*(aTemp - proposalMean)*invProposalCov*(aTemp - proposalMean)';


        %Compute densities for inverse step, needed for Metropolis
        inverseProposalMean = aTemp  + .5*optim.beta*optim.MCMC.stepWidth^2*grad_qprop;
        invProposalExponent = -.5*(a - inverseProposalMean)*invProposalCov*(a - inverseProposalMean)';
        aDiffTemp = (aTemp - conductivity.mu_a);
        pExponentTemp = -.5*optim.beta*aDiffTemp*conductivity.covar_aInv*aDiffTemp';

        Metropolis = exp(invProposalExponent - proposalExponent + ...
            optim.beta*logUprop - optim.beta*logU + pExponent - pExponentTemp);
    elseif(strcmp(optim.samplingMethod,'nonlocal'))
        %propose from p(a|theta) directly
        aTemp = mvnrnd(conductivity.mu_a, conductivity.covar_a);
        
        lambda = computeLambda(aTemp, domain, conductivity.lambdaCutoff);
        
        [logUprop] = outFunction(lambda, aTemp, conductivity, physical, domain);
        
        pExponentTemp = 0; %dummy
        
        Metropolis = exp(logUprop - logU);
    else
        error('Unknown MCMC method')
    end
    
    r = rand;
    if(r < Metropolis)
       pExponent = pExponentTemp;
       logU = logUprop;
       gradLogU = gradLogUprop;
       grad_q = grad_qprop;
       a = aTemp;
       accepted_therm = accepted_therm + 1;
    end
end
if(optim.MCMC.nTherm)
    accThermRatio = accepted_therm/optim.MCMC.nTherm;
    refineWidth = 1;
else
    %no refinement of step width
    refineWidth = 0;
end

%After thermalization, a refinement of the step width should be allowed
%Refine step width according to acceptance ratio, tune it to about .7
% accThermRatio
if(accThermRatio && refineWidth)
%     oldStepWidth = optim.MCMC.stepWidth
    optim.MCMC.stepWidth = (1/.7)*accThermRatio*optim.MCMC.stepWidth;
%     newStepWidth = optim.MCMC.stepWidth
elseif(refineWidth)
    accThermRatio
    warning('zero acceptance in thermalization, take .5 of stepWidth')
    optim.MCMC.stepWidth = .5*optim.MCMC.stepWidth;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%actual sampling
accepted = 0;
currentAcceptance = 0;
j = 1;
for i = 1:nSamplesIteration
    if(strcmp(optim.samplingMethod,'langevin'))
        %Langevin proposal 15/02/16
        proposalMean = a + .5*optim.beta*optim.MCMC.stepWidth^2*grad_q;
        %outFunctionOut.gradient_a
        aTemp = proposalMean + optim.MCMC.stepWidth*mvnrnd(zMean, zCovariance);
        lambda = computeLambda(aTemp, domain, conductivity.lambdaCutoff);

        [logUprop, gradLogUprop] = outFunction(lambda, aTemp, conductivity, physical, domain);
        grad_qprop = -(conductivity.covar_a\((aTemp - conductivity.mu_a)'))' + gradLogUprop;

        proposalExponent = -.5*(aTemp - proposalMean)*invProposalCov*(aTemp - proposalMean)';

        %Compute densities for inverse step, needed for Metropolis
        inverseProposalMean = aTemp  + .5*optim.beta*optim.MCMC.stepWidth^2*grad_qprop;
        invProposalExponent = -.5*(a - inverseProposalMean)*invProposalCov*(a - inverseProposalMean)';
        aDiffTemp = (aTemp - conductivity.mu_a);
        pExponentTemp = -.5*optim.beta*aDiffTemp*conductivity.covar_aInv*aDiffTemp';

        Metropolis = exp(invProposalExponent - proposalExponent + ...
            optim.beta*logUprop - optim.beta*logU + pExponentTemp - pExponent);

    elseif(strcmp(optim.samplingMethod,'nonlocal'))
        %propose from p(a|theta) directly
        aTemp = mvnrnd(conductivity.mu_a, conductivity.covar_a);
        
        lambda = computeLambda(aTemp, domain, conductivity.lambdaCutoff);
        
        [logUprop] = outFunction(lambda, aTemp, conductivity, physical, domain);
        
        pExponentTemp = 0; %dummy
        
        Metropolis = exp(logUprop - logU);
    else
        error('Unknown MCMC sampling method')
    end
    
    r = rand;
    if(r < Metropolis)
       pExponent = pExponentTemp;
       logU = logUprop;
       gradLogU = gradLogUprop;
       grad_q = grad_qprop;
       a = aTemp;
       accepted = accepted + 1;
       currentAcceptance = accepted/i;
    end
    
    if(mod(i,optim.MCMC.nGap) == 0)
        statistics.logUmean = (j/(j + 1))*statistics.logUmean + (1/(j + 1))*logU;
        p = mvnpdf(a,conductivity.mu_a,conductivity.covar_a);
        q = exp(logU)*p;
        statistics.q_mean = (j/(j + 1))*statistics.q_mean + (1/(j + 1))*q;
        statistics.grad_q_mean = (j/(j + 1))*statistics.grad_q_mean + (1/(j + 1))*grad_q;
        statistics.gradLogUMean = (j/(j + 1))*statistics.gradLogUMean + (1/(j + 1))*gradLogU;

        samples(j,:) = a;
        j = j + 1;
    end
    
    if(mod(i,100) == 0)
        iterationProgress = i/nSamplesIteration
        currentAcceptance
        %break if we exceed runtime limit
        if(toc > optim.MCMC.runTimeLimit)
            warning('Exceeded run time limit for sampling - stopped sampling')
            break;
        end

    end
        
end

if(j <= sampleSize)
    if(isnan(samples(j,1)))
        samples(j:end,:) = [];
    end
end
acceptance = accepted/nSamplesIteration;
% Out.acceptance_ratio
if(acceptance < .1)
    warning('Metropolis acceptance dropped below .1')
end

%wrap the output in a structure
out.samples = samples;
out.acceptance = acceptance;
out.optim = optim;
out.logU = logU;
out.pExponent = pExponent;
out.statistics = statistics;


end


