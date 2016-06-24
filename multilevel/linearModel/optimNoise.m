function [optimSigma] = optimNoise(output, designMatrix, priorMean, priorCov, mu_sigma)
%Find max marginal likelihood noise parameter sigma

% %negMargLogPost is the negative marginal logarithmic posterior, i.e. nL =
% %-log(Prod_i p(u_i|sigma)p(sigma)) (up to constants). p(sigma) = 1/mu_sigma
% %e^(-sigma/mu_sigma)
%     function nL = negMargLogPost(logSigma)
%         sigma = exp(logSigma);
%         nL = 0;
%         for i = 1:size(output,1)
%             phi = [1, input(i)];
%             nL = nL + log(sigma^2 + phi*postCov*phi') +...
%                 + ((output(i) - phi*postMean)^2)/...
%                 ((sigma^2 + phi*postCov*phi')^2) + sigma/mu_sigma;
%         end
%     end
% 
% options = optimoptions('fminunc','Algorithm','quasi-newton');
% fun = @(logSig)negMargLogPost(logSig);
% logSigmaOpt = fminunc(fun, 0,options);
% 
% optimSigma = exp(logSigmaOpt);


%negMargLogPost is the negative marginal logarithmic posterior, i.e. nL =
%-log(p(\vec u_s|sigma)) - log(p(sigma))
    function [nL, gnL] = negMargLogPost(logSigma)
        sigma = exp(logSigma);

        A = sigma^2*eye(size(output,1)) + designMatrix*priorCov*designMatrix';
        invA = inv(A);
        v = output - designMatrix*priorMean;
%         detA = det(A);
        %compute determinant via eigenvalues for safety
        logDetA = logdet(A);
        logDetA = real(logDetA);
%         nL = .5*logDetA + .5*v'*(invA*v) + sigma/mu_sigma;
        nL = .5*logDetA + .5*v'*(A\v) + sigma/mu_sigma;
%         gnL = sigma^2*trace(invA) - sigma^2*v'*invA*invA*v + sigma/mu_sigma;
        gnL = sigma^2*trace(invA) - sigma^2*(v'/A)*(A\v) + sigma/mu_sigma;
    end

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','Display','iter','MaxFunEvals',400);
fun = @(logSig)negMargLogPost(logSig);

[logSigmaOpt, val] = fminunc(fun, 0,options);



optimSigma = exp(logSigmaOpt);

end

