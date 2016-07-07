function [F] = get_glob_force_gradientGP(domain, k, e)
    %Assemble global force vector
    
    F = zeros(domain.nEquations,1);
    
    f = get_loc_force_gradientGP(domain, k(:,:,e), e);
    for ln = 1:4
        eqn = domain.lm(e,ln);
        if(eqn)
            F(eqn) = F(eqn) + f(ln);
        end
    end
    
end