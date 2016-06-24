function [out_val_mat] = basis_functions(conductivity, domain)
%Basis functions for cont lambda

%The output of the basis functions at specific coordinates. Each FEM
%element has its values in a row of out_val_mat
out_val_mat = zeros(domain.nElements,conductivity.Nx*conductivity.Ny);

%polynomials
k = 1;

%old version
% for i = 0:(conductivity.Nx - 1)
%     for j = 0:(conductivity.Ny - 1)
%         for e = 1:domain.nElements
%            out_val_mat(e,k) = (domain.coord_mat(1,e)^i)*(domain.coord_mat(2,e)^j);
%         end
%         k = k + 1;
%     end
% end

if(strcmp(domain.basisFunctionType, 'polynomial')) %polynomial basis functions
    for i = 0:(conductivity.Nx - 1)
        for j = 0:(conductivity.Ny - 1)
            for e = 1:domain.nElements
               out_val_mat(e,k) = (domain.coord_mat(1,e)^i)*(domain.coord_mat(2,e)^j);
            end
            k = k + 1;
        end
    end
elseif(strcmp(domain.basisFunctionType,'gauss')) %gaussian basis functions
    for i = 0:(conductivity.Nx - 1)
        for j = 0:(conductivity.Ny - 1)
            %xDiff = x - mu for each element
            xDiff = domain.coord_mat;
            xDiff = xDiff(1:2,:);
            gaussLocation = repmat(domain.gaussLocations(k,:)',1,domain.nElements);
            xDiff = xDiff - gaussLocation;
            exponent = (-1/(2*domain.gaussVariance))*...
                diag(xDiff'*xDiff);

            out_val_mat(:,k) = exp(exponent);
            k = k + 1;
        end
    end
else
    error('Unknown lambda basis function type')
end