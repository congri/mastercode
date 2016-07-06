function [log_q, d_log_q] = integrand(input, outFunction, conductivity, physical, domain, beta)
%Integrand to sample from, i.e. q(z) = U(z)*p(z)
    
    

[logU, gradLogU] = outFunction(input, conductivity, physical, domain);

inDiff = (input - conductivity.mu_a);
logP = -.5*inDiff*conductivity.covar_aInv*inDiff';
gradLogP = -(conductivity.covar_aInv*((input - conductivity.mu_a)'))';

log_q = beta*(logU + logP);
d_log_q = beta*(gradLogU + gradLogP);
end

