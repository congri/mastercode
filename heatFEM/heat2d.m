function [Out] = heat2d(domain, physical, control, D)
%2D heat conduction main function
%Gives back temperature on point x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nested function for performance; B-matrices should be constant for each
%element
%short hand notation
x1 = domain.lc(1,1,1);
x2 = domain.lc(1,2,1);
y1 = domain.lc(1,1,2);
y4 = domain.lc(1,4,2);

%Gauss points
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
eta1 = -1/sqrt(3);
eta2 = 1/sqrt(3);
    
%Coordinate transformation
xI = 0.5*(x1 + x2) + 0.5*xi1*(x2 - x1);
xII = 0.5*(x1 + x2) + 0.5*xi2*(x2 - x1);
yI = 0.5*(y1 + y4) + 0.5*eta1*(y4 - y1);
yII = 0.5*(y1 + y4) + 0.5*eta2*(y4 - y1);

%Assuming bilinear shape functions here!!!
B1 = (1/4)*[yI-y4 y4-yI yI-y1 y1-yI; xI-x2 x1-xI xI-x1 x2-xI];
B2 = (1/4)*[yII-y4 y4-yII yII-y1 y1-yII; xII-x2 x1-xII xII-x1 x2-xII];
%Do not forget cross terms
B3 = (1/4)*[yI-y4 y4-yI yI-y1 y1-yI; xII-x2 x1-xII xII-x1 x2-xII];
B4 = (1/4)*[yII-y4 y4-yII yII-y1 y1-yII; xI-x2 x1-xI xI-x1 x2-xI];
Bvec = [B1; B2; B3; B4];
BvecT = Bvec';

ivec = [1; 2; 1; 2; 3; 4; 3; 4; 5; 6; 5; 6; 7; 8; 7; 8];
jvec = [1; 1; 2; 2; 3; 3; 4; 4; 5; 5; 6; 6; 7; 7; 8; 8];
function [k] = get_loc_stiff2(D)
%Gives the local stiffness matrix
vvec = repmat(D(:),4,1);

Dmat = sparse(ivec, jvec, vvec, 8, 8);

k = BvecT*Dmat*Bvec;

end





%Compute local stiffness matrices, once and for all
Out.localStiffness = zeros(4,4,domain.nElements);
for e = 1:domain.nElements
    Out.localStiffness(:,:,e) = get_loc_stiff2(D(:,:,e));
end

%Global stiffness matrix
Out.globalStiffness = get_glob_stiff(domain,Out.localStiffness);
%Global force vector
Out.globalForce = get_glob_force(domain,physical,Out.localStiffness);

%Finally solving the equation system
Out.naturalTemperatures = Out.globalStiffness\Out.globalForce;

%Temperature field
Tf = zeros(domain.totalNodeNumber,1);

Tf(domain.id) = Out.naturalTemperatures;
if(strcmp(domain.lowerBoundary,'essential'))
    Tf(1:(domain.Nx + 1)) = physical.Tb(1);
end
if(strcmp(domain.rightBoundary,'essential'))
   Tf((domain.Nx + 1):(domain.Nx + 1):(domain.Nx + 1)*(domain.Ny + 1)) = physical.Tb(2);
end
if(strcmp(domain.upperBoundary,'essential'))
    Tf(domain.Ny*(domain.Nx + 1):(domain.Ny + 1)*(domain.Nx + 1)) = physical.Tb(3);
end
if(strcmp(domain.leftBoundary,'essential'))
   Tf(1:(domain.Nx + 1):(domain.Nx + 1)*(domain.Ny + 1)) = physical.Tb(4); 
end

Tff = zeros(domain.Nx + 1, domain.Ny + 1);

for i = 1:domain.totalNodeNumber
   Tff(i) = Tf(i);
   if(~isnan(domain.nodalCoordinates(4,i)))
       Tff(i) = domain.nodalCoordinates(4,i);
   end
end
Tff = Tff';
Out.Tff = Tff;

%Find temperature Tx on input point x
[X,Y] = meshgrid(linspace(0,domain.length,domain.Nx + 1), linspace(0, domain.height, domain.Ny + 1));
%'nearest' might be changed by 'linear' for better accuracy. However, the
%gradient becomes much simpler with 'nearest' interpolation
Out.temperatureMeasurements = interp2(X,Y,Tff,physical.x(:,1),physical.x(:,2),'nearest');

%Plot?
if(control.plt)
    plotHeatMap;
    %[Qx, Qy] = plotflux(X,Y,Tff,domain.elementLength,domain.Nx,domain.Ny,D);
end
end