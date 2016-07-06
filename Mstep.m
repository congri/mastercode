function [paramsOpt] = Mstep(samples, K_dK, paramsOpt)
%M-step of Expectation-Maximization
    

Nsamples = size(samples, 1);

%Robbins-Monro
[K, dK] = K_dK(paramsOpt);

Kinv_dK = 0*dK;
dL = zeros(1, size(paramsOpt, 2));


%To estimate Robbins-Monro stepsize
for i = 1:size(paramsOpt, 2)
    Kinv_dK(:,:,i) = K\dK(:,:,i);
    dL(i) = -.5*trace(Kinv_dK(:,:,i)) + (1/(2*Nsamples))*trace(samples*(Kinv_dK(:,:,i)/K)*samples');
end

RMoffset = 30;
stepSize = [RMoffset/norm(dL) RMoffset/norm(dL)];

converged = false;
iterations = 0;
maxIterations = 200;
convThreshold = 1e-4;
while(~converged)
    dLold = dL;
    for i = 1:size(paramsOpt, 2)
        Kinv_dK(:,:,i) = K\dK(:,:,i);
        dL(i) = -.5*trace(Kinv_dK(:,:,i)) + (1/(2*Nsamples))*trace(samples*(Kinv_dK(:,:,i)/K)*samples');
    end
    
%     dL
    if(norm(dL) > 1e4)
        dL = (1e4/norm(dL))*dL;
    end
    for i = 1:size(paramsOpt, 2)
        %decrease stepSize if there is a sign change
        if(dL(i)*dLold(i) < 0)
            stepSize(i) = .8*stepSize(i);
        end       
    end
    %slightly increase stepsize if gradient is not changing a lot
    if(norm(dL - dLold)/norm(dL) < 1e-3)
        dL = 1.1*dL;
    end
    dLreg = dL
    
    
    paramsOpt = paramsOpt + (stepSize/(RMoffset + iterations)).*dL;
    
    [K, dK] = K_dK(paramsOpt);
    iterations = iterations + 1;

    ndL = norm(dL)
    if(iterations > maxIterations || norm(dL) < convThreshold)
        converged = true;
    end
end
paramsOpt
dL

end

