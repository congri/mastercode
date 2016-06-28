%FEM domain specific parameters

domain.length = 1;      %length of domain
domain.height = 1;      %height of domain
%boundary conditions
%edge length of element
domain.elementLength = sqrt(domain.length*domain.height/domain.nElements);
if(mod(domain.length,domain.elementLength) ~= 0 || mod(domain.height,domain.elementLength ~= 0))
    error('Linear element number is not a divisor of edge length')
end
%number of elements in x and y direction
domain.Nx = domain.length/domain.elementLength;
domain.Ny = domain.height/domain.elementLength;
%total node number
domain.totalNodeNumber = (domain.Nx + 1)*(domain.Ny + 1);
%Nodal coordiante array
%holds global nodal coordinates in the first two lines (x and y).
%In the third line, the equation number is stored
domain.nodalCoordinates = get_coord(domain, physical);
%globalNodeNumber holds the global node number, given the element number as row
%and the local node number as column indices
domain.globalNodeNumber = get_glob(domain);
Sg = size(domain.globalNodeNumber);
%lm takes element number as row and local node number as column index
%and gives equation number
domain.lm = domain.globalNodeNumber;
for i = 1:Sg(1)
    for j = 1:Sg(2)
        domain.lm(i,j) = domain.nodalCoordinates(3,domain.globalNodeNumber(i,j));
    end
end
clear Sg i j;
%equation number and local node number precomputation for sparse stiffness
%assembly
[domain.Equations, domain.LocalNode] = get_equations(domain.nElements, domain.lm);
domain.kIndex = sub2ind([4 4 domain.nElements], domain.LocalNode(:,1), domain.LocalNode(:,2), domain.LocalNode(:,3));
domain.essNode = false(domain.nElements, 4);
for e = 1:domain.nElements
    for i = 1:4
        domain.essNode(e, i) = ~isnan(domain.nodalCoordinates(4,domain.globalNodeNumber(e,i)));
    end
end
%lc gives node coordinates, taking in element number and local node
%number
domain.lc = get_loc_coord(domain);
%Get mapping from equation number back to global node number
domain.id = get_id(domain.nodalCoordinates);
%Coordinates of centers of finite elements
domain.coord_mat = get_el_coord(domain);
domain.coordMatSquared = domain.coord_mat(1,:)'*domain.coord_mat(2,:);
domain.nEquations = max(domain.nodalCoordinates(3,:));
%Total number of equations
neq = double(max(domain.nodalCoordinates(3,:)));
domain.Kinit = spalloc(neq,neq,3*neq);
for k = 1:domain.nElements
    for i = 1:4
        %essential boundary (yes or no) given local node and element number
        domain.essentialBoundary(i,k) = ~isnan(domain.nodalCoordinates(4,domain.globalNodeNumber(k,i)));
    end
end
clear k i neq;

domain.basisFunctionType = 'gauss';    %polynomial or gauss