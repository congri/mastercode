function [F] = get_glob_force(domain, physical, k)
%Assemble global force vector

neq = max(domain.nodalCoordinates(3,:));
F = zeros(neq,1);

for i = 1:domain.nElements
    f = get_loc_force(i, domain, k, physical);
    for ln = 1:4
       eqn = domain.lm(i,ln);
       if(eqn ~= 0)
       F(eqn) = F(eqn) + f(ln);
       end
    end
end

end

