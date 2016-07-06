function [lambda, reshapedLambda] = GPconductivity(m, K)
%Draws a random heat conductivity field from a Gaussian process prior

y = mvnrnd(m, K);
lambda = exp(y);
reshapedLambda = reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda)));

end

