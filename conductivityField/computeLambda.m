function [lambda] = computeLambda(z, domain, cutoff)
%Compute heat conductivities lambda from a

if(strcmp(domain.basisFunctionType, 'polynomial'))
    lambda = exp(-(z*domain.basisFunctionValues'))' + 1e-5;
    %cutoff of lambda
    lambda(lambda > cutoff) = cutoff;
elseif(strcmp(domain.basisFunctionType, 'gauss'))
    lambda = (exp(z)*domain.basisFunctionValues')' + 1e-5;
    %cutoff of lambda
    lambda(lambda > cutoff) = cutoff;
elseif(strcmp(domain.basisFunctionType, 'GP'))
    lambda = exp(z') + 1e-5;
    %cutoff of lambda
    lambda(lambda > cutoff) = cutoff;
else
    error('unknown basis function type for lambda')
end


end

