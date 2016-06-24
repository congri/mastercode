function [predMean, predVar] = predParams(postMean, postCov, input, sigmaNoise, bFunH)
%Predictive distribution of Gaussian linear model

% predMean = zeros(size(input,1),1);
predVar = zeros(size(input,1),1);
A = bFunH(input)';
predMean = A*postMean;
for i = 1:size(input,1)
    predVar(i) = sigmaNoise^2 + A(i,:)*postCov*A(i,:)';
end

end

