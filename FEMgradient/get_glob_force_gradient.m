function [F] = get_glob_force_gradient(domain, k)
%Assemble global force vector

F = zeros(domain.nEquations,1);

f = get_loc_force_gradient(domain, k);
for e = 1:domain.nElements
    for ln = 1:4
       eqn = domain.lm(e,ln);
       if(eqn)
          F(eqn) = F(eqn) + f(ln,e);
       end
    end
end

end