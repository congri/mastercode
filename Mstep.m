function [paramsOpt] = Mstep(samples, K_dK, paramsOpt)
    %M-step of Expectation-Maximization
    
    
    Nsamples = size(samples, 1);
    
    %Robbins-Monro
    [K, grad] = K_dK(paramsOpt);
    
    dL = zeros(1, size(paramsOpt, 2));
    
    
    Kinv_dK = grad; %preallocation
    %To estimate Robbins-Monro stepsize
    for i = 1:size(paramsOpt, 2)
        Kinv_dK(i).KidK = K\grad(i).dK;
        dL(i) = -.5*trace(Kinv_dK(i).KidK) + (1/(2*Nsamples))*trace(samples*(Kinv_dK(i).KidK/K)*samples');
    end
    
    RMoffset = 30;
    if norm(dL) > 0
        stepSize = .1*[RMoffset/norm(dL) RMoffset/norm(dL)];
    else
        stepSize = [1 1];
    end
    
    converged = false;
    iterations = 0;
    maxIterations = 200;
    convThreshold = 1e-2;
    while(~converged)
        dLold = dL;
        for i = 1:size(paramsOpt, 2)
            Kinv_dK(i).KidK = K\grad(i).dK;
            dL(i) = -.5*trace(Kinv_dK(i).KidK) + (1/(2*Nsamples))*trace(samples*(Kinv_dK(i).KidK/K)*samples');
        end
        
        for i = 1:size(paramsOpt, 2)
            %decrease stepSize if there is a sign change
            if(dL(i)*dLold(i) < 0)
                stepSize(i) = .5*stepSize(i);
            end
        end
        %slightly increase stepsize if gradient is not changing a lot
%         gradChange = norm(dL - dLold)/norm(dL)
        if(norm(dL - dLold)/norm(dL) < 1.5)
            stepSize = 1.3*stepSize;
        end
%         stepSize
        dL

        
        
        step = (stepSize/(RMoffset + iterations)).*dL;
        %Cutoff in stepsize
        cutoff = 1;
        if norm(step) > cutoff
            step = (cutoff/norm(step))*step;
        end
        paramsOpt = paramsOpt + step;
        
        [K, grad] = K_dK(paramsOpt);
        iterations = iterations + 1;
        
%         ndL = norm(dL)
        if(iterations > maxIterations || norm(dL) < convThreshold)
            converged = true;
        end
    end
    paramsOpt
    dL
    
end

