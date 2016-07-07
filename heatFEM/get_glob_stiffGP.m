function [K] = get_glob_stiffGP(domain, k, e)
%Gives global stiffness matrix K gradient only for Gaussian process
%conductivity field version!

%Can this be done more efficiently?
%Assign global stiffness matrix
% localNodeInit = 1:4;
% equations = domain.lm(l,localNodeInit);
% equations = equations(equations > 0);
% localNode = localNodeInit(equations > 0);
% [eq1, eq2] = meshgrid(equations);
K = sparse(domain.GPequations(e).eq1, domain.GPequations(e).eq2,...
 k(domain.GPequations(e).localNode,domain.GPequations(e).localNode),...
 domain.nEquations, domain.nEquations, 16);


end

