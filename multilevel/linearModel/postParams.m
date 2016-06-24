function [postMean, postCovInv] = postParams(priorMean, priorCovInv, sigmaNoise, designMatrix, dataOut)
%Compute posterior mean and covariance for a Gaussian linear model

postCovInv = priorCovInv + (sigmaNoise^(-2))*(designMatrix'*designMatrix);
postMean = postCovInv\(priorCovInv*priorMean + (sigmaNoise^(-2))*designMatrix'*dataOut);

end

