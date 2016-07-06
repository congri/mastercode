%has to be executed after conductivityParams!

%domain 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basis function values at FEM coordinate values
%for gaussian  basis functions
xG = linspace(1/(2*sqrt(conductivity.dim)), 1 - 1/(2*sqrt(conductivity.dim)), sqrt(conductivity.dim))';
xG = repmat(xG,sqrt(conductivity.dim),1);
yG = (ceil((1:conductivity.dim)/sqrt(conductivity.dim))/sqrt(conductivity.dim) - 1/(2*sqrt(conductivity.dim)))';
domain.gaussLocations = [xG yG];
domain.gaussVariance = ((1/4)*(1/conductivity.Nx))^2; %ignore if basis function type is polynomial
if(~strcmp(domain.basisFunctionType, 'GP'))
    domain.basisFunctionValues = basis_functions(conductivity, domain);
end
clear xG yG;

xDiffSq = zeros(1,domain.nElements);
for l = 1:conductivity.dim
    %xDiff = x - mu for each element
    xDiff = domain.coord_mat;
    xDiff = xDiff(1:2,:);
    gaussLocation = repmat(domain.gaussLocations(l,:)',1,domain.nElements);
    xDiff = xDiff - gaussLocation;
    for e = 1:domain.nElements
       xDiffSq(e) = xDiff(:,e)'*xDiff(:,e); 
    end
    domain.exponent(:,l) = ((-1/(2*domain.gaussVariance))*xDiffSq)';
end
clear xDiff xDiffSq l e gaussLocation;