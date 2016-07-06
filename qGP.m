function [log_q, d_log_q] = qGP(input, outFunction, conductivity, physical, domain, beta, Kinv)
%Integrand to sample from, i.e. q(z) = U(z)*p(z) 

[logU, gradLogU] = outFunction(input, conductivity, physical, domain);

logP = -.5*input*Kinv*input';
gradLogP = -input*Kinv;

log_q = beta*(logU + logP);
d_log_q = beta*(gradLogU + gradLogP);
end