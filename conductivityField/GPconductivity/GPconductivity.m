function [lambda, reshapedLambda] = GPconductivity(mu, k, x)
%Draws a random heat conductivity field from a Gaussian process prior
    

%mean vector
m = mu(x);
%covariance matrix
K = k(x);

y = mvnrnd(m, K);
lambda = exp(y);
reshapedLambda = reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda)));

end

