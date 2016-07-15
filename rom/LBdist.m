function [log_q_i] = LBdist(X, y_i, Phi, theta_c, theta_cf, physical, domain, conductivity)
    %Unnormalized lower bound distribution q_i in EM in reduced order model; see
    %Bayesian_CG_PDE.pdf
    %
    %Input:     X:          coarse model conductivity input
    %           y_i:        fine model output data point
    %           Phi:        Design matrix in p_c
    %           theta_c:    params of p_c
    %           theta_cf:   params in p_cf
    %
    %Output:    log_q_i:    log q_i = log p_cf + log p_c
    
    %Values of heat conductivities on finite elements
    [Y] = romOutput(X, conductivity, physical, domain);
    
%     size(y_i)
%     theta_cf
%     size(Y)
%     theta_c
%     size(X)
%     size(Phi)
%     size(Phi'*theta_c.theta)
%     size(X)
    
    log_q_i = -.5*(y_i - theta_cf.mu - theta_cf.W*Y)'*theta_cf.Sinv*...
        (y_i - theta_cf.mu - theta_cf.W*Y) - .5*theta_c.varInv*...
        norm(X - theta_c.theta'*Phi)^2;
    
end

