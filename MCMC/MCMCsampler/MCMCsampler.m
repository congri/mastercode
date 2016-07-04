function [out] = MCMCsampler(log_distribution, startValue, opts)
%Standalone MCMC sampler. Pass distribution, starting value and options and
%get samples
%By C. Grigo, July 2016

% log_distribution:         function handle to log probability distribution and its
%                           gradient
% startValue:               Initial value of Markov chain
% opts:                     MCMC sampling options structure

% opts:
% opts.method               proposal type: randomWalk, nonlocal or MALA
% opts.nThermalization      thermalization steps
% opts.nSamples             number of samples
% 
% %only for random walk
% 
% %only for nonlocal
% opts.nonlocal.rnd         random number generator for nonlocal proposals
% opts.nonlocal.pdf         corresponding proposal pdf
% opts.nonlocal.propMean    proposal pdf is Gaussian: mean is propMean
% opts.nonlocal.propCov     proposal pdf is Gaussian: cov is propCov

% only for MALA
% opts.MALA.stepWidth       step size parameter


rng('shuffle');     %random number seed based on system time


%preallocation of samples array
out.samples = zeros(opts.nSamples, size(startValue, 2));
samplesTherm = zeros(opts.nThermalization, size(startValue, 2));
samplesTherm(1, :) = startValue;

x = startValue;

accepted = 0;

if(strcmp(opts.method, 'MALA'))
    
   zeroMean = zeros(1, size(x, 2));
   unitCov = eye(size(x, 2));
   invProposalCov = inv(opts.MALA.stepWidth^2*unitCov);
   [log_p, d_log_p] = log_distribution(x);

else
    
    log_p = log_distribution(x);

end


%Thermalization
for i = 1:(opts.nThermalization - 1)

    if(strcmp(opts.method, 'randomWalk'))
        %Gaussian random walk MCMC
        
        xProp = mvnrnd(x, opts.randomWalk.proposalCov);
        log_pProp = log_distribution(xProp);
        Metropolis = exp(log_pProp - log_p);
        
    elseif(strcmp(opts.method, 'nonlocal'))
        %"Nonlocal" proposal distribution
        
        xProp = opts.nonlocal.rnd(opts.nonlocal.propMean, opts.nonlocal.propCov);
        log_pProp = log_distribution(xProp);
        Metropolis = exp(log_pProp - log_p)*...
            ((opts.nonlocal.pdf(x, opts.nonlocal.propMean, opts.nonlocal.propCov))/...
            (opts.nonlocal.pdf(xProp, opts.nonlocal.propMean, opts.nonlocal.propCov)));
        
    elseif(strcmp(opts.method, 'MALA'))
        %Metropolis adjusted Langevin algorithm
        
        proposalMean = x + .5*opts.MALA.stepWidth^2*d_log_p;
        xProp = proposalMean + opts.MALA.stepWidth*mvnrnd(zeroMean, unitCov);
        proposalExponent = -.5*(xProp - proposalMean)*invProposalCov*(xProp - proposalMean)';
        
        [log_pProp, d_log_pProp] = log_distribution(xProp);
        inverseProposalMean = xProp  + .5*opts.MALA.stepWidth^2*d_log_pProp;
        invProposalExponent = -.5*(x - inverseProposalMean)*invProposalCov*(x - inverseProposalMean)';
        
        Metropolis = exp(invProposalExponent - proposalExponent + log_pProp - log_p);
     
    else
        
        error('unknown MCMC sampling method')
        
    end

    r = rand;
    if(r < Metropolis)
        %Metropolis acceptance. Go to xProp

        x = xProp;
        log_p = log_pProp;
        if(strcmp(opts.method, 'MALA'))
            d_log_p = d_log_pProp;
        end
        accepted = accepted + 1;

    end
    samplesTherm(i + 1, :) = x;

end

acceptance = accepted/opts.nThermalization;
%refine proposal params after thermalization
if(strcmp(opts.method, 'randomWalk'))
    
    opts.randomWalk.proposalCov = (1/.7)*acceptance*opts.randomWalk.proposalCov;
    
elseif(strcmp(opts.method, 'nonlocal'))
    
    opts.nonlocal.propMean = mean(samplesTherm);
    %ATTENTION: fully decorrelated thermalization required!
    opts.nonlocal.propCov = .2*cov(samplesTherm) + .8*opts.nonlocal.propCov;
    
elseif(strcmp(opts.method, 'MALA'))
    %Metropolis adjusted Langevin algorithm
    
    opts.MALA.stepWidth = (1/.7)*acceptance*opts.MALA.stepWidth;

else
    error('unknown MCMC sampling method')
end



%Actual sampling
accepted = 0;
for i = 1:opts.nSamples

    if(strcmp(opts.method, 'randomWalk'))
        %Gaussian random walk MCMC

        xProp = mvnrnd(x, opts.randomWalk.proposalCov);
        Metropolis = exp(log_distribution(xProp) - log_distribution(x));

    elseif(strcmp(opts.method, 'nonlocal'))
        %"Nonlocal" proposal distribution
        
        xProp = opts.nonlocal.rnd(opts.nonlocal.propMean, opts.nonlocal.propCov);
        log_pProp = log_distribution(xProp);
        Metropolis = exp(log_pProp - log_p)*...
            ((opts.nonlocal.pdf(x, opts.nonlocal.propMean, opts.nonlocal.propCov))/...
            (opts.nonlocal.pdf(xProp, opts.nonlocal.propMean, opts.nonlocal.propCov)));
        
        
    elseif(strcmp(opts.method, 'MALA'))
        %Metropolis adjusted Langevin algorithm
        
        proposalMean = x + .5*opts.MALA.stepWidth^2*d_log_p;
        xProp = proposalMean + opts.MALA.stepWidth*mvnrnd(zeroMean, unitCov);
        proposalExponent = -.5*(xProp - proposalMean)*invProposalCov*(xProp - proposalMean)';
        
        [log_pProp, d_log_pProp] = log_distribution(xProp);
        inverseProposalMean = xProp  + .5*opts.MALA.stepWidth^2*d_log_pProp;
        invProposalExponent = -.5*(x - inverseProposalMean)*invProposalCov*(x - inverseProposalMean)';
        
        Metropolis = exp(invProposalExponent - proposalExponent + log_pProp - log_p);
        
    else
        
        error('unknown MCMC sampling method')
        
    end

    r = rand;
    if(r < Metropolis)
        %Metropolis acceptance. Go to xProp

        x = xProp;
        log_p = log_pProp;
        if(strcmp(opts.method, 'MALA'))
            d_log_p = d_log_pProp;
        end
        accepted = accepted + 1;

    end
    out.samples(i, :) = x;

end

out.acceptance = accepted/opts.nSamples;




end

