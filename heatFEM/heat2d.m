function [Out] = heat2d(domain, physical, control, D)
%2D heat conduction main function
%Gives back temperature on point x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute local stiffness matrices, once and for all
Out.localStiffness = zeros(4,4,domain.nElements);
for i = 1:domain.nElements
    Out.localStiffness(:,:,i) = get_loc_stiff(i,domain.lc,D);
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