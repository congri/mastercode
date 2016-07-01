function [p, d_log_p] = testDistribution(x, m, S)
%A probability distribution function and the gradient of its log for testing

p = mvnpdf(x, m, S);

d_log_p = - S\(x - repmat(m, size(x,1), 1))';




    
end

