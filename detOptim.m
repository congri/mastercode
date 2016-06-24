function [optim_a] = detOptim(conductivity, physical, domain)

    function [log_U, dlogU_da] = outFun(a)
        [lambda] = computeLambda(a, domain, conductivity.lambdaCutoff);
        [Out] = cont_ref_output(lambda, a, conductivity, physical, domain);
        %minus signs due to minimization/maximization
        dlogU_da = -Out.dlogU_da;
        log_U = -Out.exponent;
    end

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','Display','iter','MaxIter',100);
a0 = mvnrnd(zeros(1,conductivity.dim),eye(conductivity.dim));

fun = @outFun;
optim_a = fminunc(fun,a0,options);


end

