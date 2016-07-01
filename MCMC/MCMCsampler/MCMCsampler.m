function [samples, acceptanceRatio] = MCMCsampler(distribution, startValue, opts)
%Standalone MCMC sampler. Pass distribution, starting value and options and
%get samples
%By C. Grigo, July 2016

%   distribution:   function handle to probability distribution and its
%                   gradient
%   startValue:     Initial value of Markov chain
%   opts:           MCMC sampling options structure


rng('shuffle');     %random number seed based on system time


%preallocation of samples array
samples = zeros(opts.nSamples, size(startValue, 2));
samples(1, :) = startValue;
x = startValue;

%!!!!!!!!!!!!!!! introduce log scale!!!!!!!!!!!!!!

if(strcmp(opts.method, 'randomWalk'))
    %random walk MCMC sampling
    
    %Thermalization
    accepted = 0;
    for i = 1:(opts.nThermalization - 1)
        
        xProp = mvnrnd(x, opts.proposalCov);
        Metropolis = distribution(xProp)/distribution(x);
        
        r = rand;
        if(r < Metropolis)
           %Metropolis acceptance. Go to xProp
           
           x = xProp;
           accepted = accepted + 1;

        end
        
    end
    
    acceptanceRatio = accepted/opts.nThermalization;
    %refine step width after thermalization
    opts.proposalCov = (1/.7)*acceptanceRatio*opts.proposalCov;
    
    
    
    accepted = 0;
    for i = 1:(opts.nSamples - 1)
        
        xProp = mvnrnd(x, opts.proposalCov);
        Metropolis = distribution(xProp)/distribution(x);
        
        r = rand;
        if(r < Metropolis)
           %Metropolis acceptance. Go to xProp
           
           x = xProp;
           accepted = accepted + 1;

        end
        samples(i + 1, :) = x;
        
    end
    
    acceptanceRatio = accepted/opts.nSamples;
    
elseif(strcmp(opts.method, 'nonlocal'))
    %nonlocal Gaussian proposal

    %Thermalization
    accepted = 0;
    for i = 1:(opts.nThermalization - 1)
        
        xProp = opts.nonlocalrnd();
        Metropolis = (distribution(xProp)*opts.nonlocalpdf(x))/(distribution(x)*opts.nonlocalpdf(xProp));
        
        r = rand;
        if(r < Metropolis)
           %Metropolis acceptance. Go to xProp
           
           x = xProp;
           accepted = accepted + 1;

        end
        
    end
    
    acceptanceRatio = accepted/opts.nThermalization;
    %refine step width after thermalization
    opts.proposalCov = (1/.7)*acceptanceRatio*opts.proposalCov;
    
    
    
    accepted = 0;
    for i = 1:(opts.nSamples - 1)
        
        xProp = mvnrnd(x, opts.proposalCov);
        Metropolis = distribution(xProp)/distribution(x);
        
        r = rand;
        if(r < Metropolis)
           %Metropolis acceptance. Go to xProp
           
           x = xProp;
           accepted = accepted + 1;

        end
        samples(i + 1, :) = x;
        
    end
    
    acceptanceRatio = accepted/opts.nSamples;
    
else
   error('unknown MCMC sampling method') 
end
    
    
end

