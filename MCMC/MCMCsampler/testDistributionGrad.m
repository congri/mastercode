function [d_log_p] = testDistributionGrad(x, m, S)
%A probability distribution function and the gradient of its log for testing



M = repmat(m, size(x,1), 1);
d_log_p = - (S\(x - M)')';

    
end
