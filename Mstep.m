function [paramsOpt] = Mstep(samples, K_dK, paramsOpt)
    %M-step of Expectation-Maximization
    
    
    Nsamples = size(samples, 1);
    
    function [gL] = grL(samp, pOpt)
        
        [C, gr] = K_dK(pOpt);
        Ki_dK = gr; %preallocation
        gL = zeros(1, size(pOpt, 2));
        for j = 1:size(pOpt, 2)
            Ki_dK(j).KidK = C\gr(j).dK;
            gL(j) = -.5*trace(Ki_dK(j).KidK) + (1/(2*Nsamples))*trace(samp*(Ki_dK(j).KidK/C)*samp');
        end
        
    end
    
    %method RM uses simple Robbins-Monro optimization code below
    %method fsolve uses Matlab built-in fsolve function (preferable)
    method = 'fsolve';
    convThreshold = 1e-2;
    
    if(strcmp(method, 'fsolve'))
        
        tic
        fun = @(pOp) grL(samples, pOp);
        options = optimoptions('fsolve','Display','iter', 'TolX', 1e-2);
        paramsOpt  = fsolve(fun, paramsOpt, options)
        t1 = toc
        
    elseif(strcmp(method, 'RM'))
        
        dL = grL(samples, paramsOpt);
        tic
        RMoffset = 30;
        if norm(dL) > 0
            stepSize = .1*[RMoffset/norm(dL) RMoffset/norm(dL)];
        else
            stepSize = [1 1];
        end
        
        converged = false;
        iterations = 0;
        maxIterations = 200;
        while(~converged)
            dLold = dL;
            dL = grL(samples, paramsOpt)
            
            for i = 1:size(paramsOpt, 2)
                %decrease stepSize if there is a sign change
                if(dL(i)*dLold(i) < 0)
                    stepSize(i) = .5*stepSize(i);
                end
            end
            %slightly increase stepsize if gradient is not changing a lot
            if(norm(dL - dLold)/norm(dL) < 1.5)
                stepSize = 1.3*stepSize;
            end
            
            step = (stepSize/(RMoffset + iterations)).*dL;
            %Cutoff in stepsize
            cutoff = 1;
            if norm(step) > cutoff
                step = (cutoff/norm(step))*step;
            end
            paramsOpt = paramsOpt + step;
            
            iterations = iterations + 1;
            
            %         ndL = norm(dL)
            if(iterations > maxIterations || norm(dL) < convThreshold)
                converged = true;
            end
        end
        paramsOpt
        dL
        t2 = toc
        
    else
        error('Unknown optimization method')
    end
    
end

