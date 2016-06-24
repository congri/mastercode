function [gradK, gradF] = gradKgradF(FEMout, conductivity, domain, a)
%Computes gradient of global stiffness matrix K and global force vector F
%w.r.t. a

grad_kl = zeros(4, 4, domain.nElements, conductivity.dim);
SK = size(FEMout.globalStiffness);
gradK = zeros(SK(1),SK(2),conductivity.dim);
gradF = zeros(SK(1), conductivity.dim);
l = 1;

if(strcmp(domain.basisFunctionType,'polynomial'))
    for i = 1:conductivity.Nx
        for j = 1:conductivity.Ny
            for k = 1:domain.nElements
                grad_kl(:,:,k,l) = - domain.coord_mat(1,k)^(i - 1)*domain.coord_mat(2,k)^(j - 1)*FEMout.localStiffness(:,:,k)...
                    + domain.coord_mat(1,k)^(i - 1)*domain.coord_mat(2,k)^(j - 1)*conductivity.localStiffnessOffset(:,:,k);
            end
            %Assemble gradient of global stiffness matrix
            gradK(:,:,l) = get_glob_stiff(domain, grad_kl(:,:,:,l));
            %Assemble gradient of global force vector
            gradF(:,l) = get_glob_force_gradient(domain, grad_kl(:,:,:,l));
            l = l + 1;
        end   
    end
elseif(strcmp(domain.basisFunctionType,'gauss'))%performance optimized
    localNodeInit = 1:4;
    for k = 1:domain.nElements
        equations = domain.lm(k,localNodeInit);
        equations = equations(equations ~= 0);
        localNode = localNodeInit(equations ~= 0);
        %Boundary value temperature of element e
        Tbflag = 0;
        Tb = zeros(4,1);
        for i = 1:4
            if(domain.essentialBoundary(i,k))
                Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(k,i));
                Tbflag = 1;
            end
        end
        prefactor = exp(a + domain.exponent(k,:)');
        for l = 1:conductivity.dim
            %stiffness matrix
            %localStiffnessOffset is with unity heat conductivity here
            grad_kl0 = prefactor(l)*conductivity.localStiffnessOffset;
            gradK(equations, equations,l) = gradK(equations, equations,l) + grad_kl0(localNode,localNode);
        
            %force vector
            if(Tbflag)
                f = -grad_kl0*Tb;
                for ln = 1:4
                   eqn = domain.lm(k,ln);
                   if(eqn)
                      gradF(eqn,l) = gradF(eqn,l) + f(ln);
                   end
                end
            end
        end
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

