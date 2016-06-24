%Conductivity field params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


conductivity.Nx = 3; %number of basis functions in x-direction
conductivity.Ny = 3; %number of basis functions in y-direction
component = 6; %optimization component

%start parameters for distribution on a
conductivity.covar_a = (1/100)*eye(conductivity.Nx*conductivity.Ny);
conductivity.covar_aInv = inv(conductivity.covar_a);
% conductivity.mu_a = [-2 -1 1 0 -1 -1 -1 0 0];
conductivity.dim = conductivity.Nx*conductivity.Ny;
conductivity.mu_a = zeros(1,conductivity.dim);
conductivity.lambdaCutoff = 1e8;
conductivity.muRegularization = 0e0;  %regularization param to avoid large mu_a; equvivalent to a prior
                                %on mu_a with precision = 2*optim.muReuglarization
%We take an offset of 1e-5 in lambda to avoid singular conductivity matrices.
%Therefore, we need an offset local conductivity matrix
D = zeros(2,2,domain.nElements);
if(strcmp(domain.basisFunctionType, 'polynomial'))
    for i = 1:domain.nElements
        D(:,:,i) = 1e-5*eye(2); %only isotropic material
    end
    for e = 1:domain.nElements
        conductivity.localStiffnessOffset(:,:,e) = get_loc_stiff(e,domain.lc,D);
    end
elseif(strcmp(domain.basisFunctionType, 'gauss'))
    for i = 1:domain.nElements
        D(:,:,i) = eye(2); %only isotropic material
    end
    conductivity.localStiffnessOffset = get_loc_stiff(1,domain.lc,D);
else
    error('Unknown basis function for lambda')
end
clear D i e;