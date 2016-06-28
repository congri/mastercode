function [K] = get_glob_stiff2(domain, k)
%Gives global stiffness matrix K

% n = size(domain.LocalNode, 1);
% kvec = zeros(n, 1);
% for i = 1:n
%     kvec(i) = k(domain.LocalNode(i, 1), domain.LocalNode(i, 2), domain.LocalNode(i,3));
% end

K = sparse(domain.Equations(:,1), domain.Equations(:,2), k(domain.kIndex));

end

