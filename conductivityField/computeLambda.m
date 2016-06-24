function [lambda] = computeLambda(a, domain, cutoff)
%Compute heat conductivities lambda from a

if(strcmp(domain.basisFunctionType, 'polynomial'))
    lambda = exp(-(a*domain.basisFunctionValues'))' + 1e-5;
    %cutoff of lambda
    lambda(lambda > cutoff) = cutoff;
elseif(strcmp(domain.basisFunctionType, 'gauss'))
    lambda = (exp(a)*domain.basisFunctionValues')' + 1e-5;
    %cutoff of lambda
    lambda(lambda > cutoff) = cutoff;
else
    error('unknown basis function type for lambda')
end


end

