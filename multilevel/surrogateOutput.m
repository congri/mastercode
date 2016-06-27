function [predExponent, grad_log_UPred, output] = surrogateOutput(lambda, input, conductivity, physical, domain,...
     eta, sigmaNoise, bFunH)
%surrogate output and gradient

if(nargout > 1)
    [logU, gradLogU] = cont_ref_output(lambda, input, conductivity, physical, domain);
else
    logU = cont_ref_output(lambda, input, conductivity, physical, domain);
end

uf = logphiInv(logU, 0);
[phi, dphi] = bFunH(uf);
uPred = eta*phi;
% Q = uPred/(sqrt(1 + sigmaNoise^2));

if(nargout > 2)
    output = normcdf(uPred);
end

predExponent = logphi(uPred);       %median prediction
% predExponent = logphi(Q);         %mean prediction

%mean prediction gradient
% grad_log_UPred = ((eta*dphi)/sqrt(1 + sigmaNoise^2))*...
%     exp(-.5*Q^2 + .5*uf^2 + logU - predExponent)*gradLogU;

%median prediction gradient
if(nargout > 1)
    grad_log_UPred = (eta*dphi)*...
        exp(-.5*uPred^2 + .5*uf^2 + logU - predExponent)*gradLogU;
end

end



