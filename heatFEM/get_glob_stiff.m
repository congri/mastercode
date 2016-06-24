function [K] = get_glob_stiff(domain,k)
%Gives the global stiffness matrix

K = domain.Kinit;

%Can this be done more efficiently?
%Assign global stiffness matrix
localNodeInit = 1:4;
for e = 1:domain.nElements
    equations = domain.lm(e,localNodeInit);
    equations = equations(equations > 0);
    localNode = localNodeInit(equations > 0);
    K(equations, equations) = K(equations, equations) + k(localNode,localNode,e);
end


end

