function [grad, gradF] = gradKgradF(FEMout, conductivity, domain, a)
%Computes gradient of global stiffness matrix K and global force vector F
%w.r.t. a

SK = size(FEMout.globalStiffness,1);
grad.K = sparse(SK,SK,3*SK);
grad = repmat(grad, 1, conductivity.dim);
gradF = zeros(SK(1), conductivity.dim);
l = 1;

if(strcmp(domain.basisFunctionType,'polynomial'))
    grad_kl = zeros(4, 4, domain.nElements, conductivity.dim);
    for i = 1:conductivity.Nx
        for j = 1:conductivity.Ny
            for k = 1:domain.nElements
                grad_kl(:,:,k,l) = - domain.coord_mat(1,k)^(i - 1)*domain.coord_mat(2,k)^(j - 1)*FEMout.localStiffness(:,:,k)...
                    + domain.coord_mat(1,k)^(i - 1)*domain.coord_mat(2,k)^(j - 1)*conductivity.localStiffnessOffset(:,:,k);
            end
            %Assemble gradient of global stiffness matrix
            grad(l).K = get_glob_stiff2(domain, grad_kl(:,:,:,l));
            %Assemble gradient of global force vector
            gradF(:,l) = get_glob_force_gradient(domain, grad_kl(:,:,:,l));
            l = l + 1;
        end   
    end
elseif(strcmp(domain.basisFunctionType,'gauss'))%performance optimized
    grad_kl = zeros(4, 4, domain.nElements);

    for l = 1:conductivity.dim
        for k = 1:domain.nElements
            grad_kl(:,:,k) = exp(a(l) + domain.exponent(k,l)')*conductivity.localStiffnessOffset;
        end
        grad(l).K = get_glob_stiff2(domain, grad_kl);
        gradF(:,l) = get_glob_force_gradient(domain, grad_kl);
    end

else
    error('unknown basis function type for lambda')
end





%old version
% for i = 1:conductivity.Nx
%     for j = 1:conductivity.Ny
%         for k = 1:domain.nElements
%             grad_kl(:,:,k,l) = - domain.coord_mat(1,k)^(i - 1)*domain.coord_mat(2,k)^(j - 1)*FEMout.localStiffness(:,:,k)...
%                 + domain.coord_mat(1,k)^(i - 1)*domain.coord_mat(2,k)^(j - 1)*conductivity.localStiffnessOffset(:,:,k);
%         end
%         %Assemble gradient of global stiffness matrix
%         gradK(:,:,l) = get_glob_stiff(domain, grad_kl(:,:,:,l));
%         %Assemble gradient of global force vector
%         gradF(:,l) = get_glob_force_gradient(domain, grad_kl(:,:,:,l));
%         l = l + 1;
%     end   
% end




end

