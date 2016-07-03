function [log_p, d_log_p] = testDistribution(x, m, S)
%A probability distribution function and the gradient of its log for testing



M = repmat(m, size(x,1), 1);
log_p = - .5*size(m, 2)*log(2*pi) - .5*logdet(S) - .5*diag((x - M)*(S\(x - M)'));

d_log_p = - S\(x - M)';

    
end

