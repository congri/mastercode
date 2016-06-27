function [nc] = get_coord(domain, physical)
%Gives nodal coordinates in the first two rows and equation number from 
%global node number in the third row. Temperature of essential boundaries
%is given in the fourth row, heat flux on natural boundaries in the fifth

%comparing precision
p = 10^(-15);

%Assign global node coordinates and equation numbers
x = 0;
y = 0;
j = 1;
nc = NaN*zeros(5,domain.totalNodeNumber);
for i = 1:domain.totalNodeNumber
    nc(1,i) = x;
    nc(2,i) = y;

%Assign equation numbers
boundaryflag = false; %false, if the node is an internal node

if(i == 1) %lower left corner
    boundaryflag = true;
    if(strcmp(domain.lowerBoundary,'essential') || strcmp(domain.leftBoundary,'essential'))
    %if any of the two boundaries is essential, the node is essential
    nc(3,i) = 0;
    if(strcmp(domain.lowerBoundary,'essential') && strcmp(domain.leftBoundary,'essential'))
    nc(4,i) = 0.5*(physical.Tb(1) + physical.Tb(4)); %We take the mean value of boundaries
    elseif(strcmp(domain.lowerBoundary,'essential'))
    nc(4,i) = physical.Tb(1);
    elseif(strcmp(domain.leftBoundary,'essential'))
    nc(4,i) = physical.Tb(4);
    else
    error('No essential boundary?');
    end
    elseif(strcmp(domain.lowerBoundary,'natural') && strcmp(domain.leftBoundary, 'natural'))
    %if both boundaries are natural, the node is natural
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = 0.5*(physical.qb(1) + physical.qb(4)); %we take the mean value
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(i == domain.Nx + 1) %lower right corner
    boundaryflag = true;
    if(strcmp(domain.lowerBoundary,'essential') || strcmp(domain.rightBoundary,'essential'))
    %if any of the two boundaries is essential, the node is essential
    nc(3,i) = 0;
    if(strcmp(domain.lowerBoundary,'essential') && strcmp(domain.rightBoundary,'essential'))
    nc(4,i) = 0.5*(physical.Tb(1) + physical.Tb(2)); %We take the mean value of boundaries
    elseif(strcmp(domain.lowerBoundary,'essential'))
    nc(4,i) = physical.Tb(1);
    elseif(strcmp(domain.rightBoundary,'essential'))
    nc(4,i) = physical.Tb(2);
    else
    error('No essential boundary?');
    end
    elseif(strcmp(domain.lowerBoundary,'natural') && strcmp(domain.rightBoundary, 'natural'))
    %if both boundaries are natural, the node is natural
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = 0.5*(physical.qb(1) + physical.qb(2)); %we take the mean value
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(i == (domain.Nx + 1)*(domain.Ny + 1)) %upper right corner
    boundaryflag = true;
    if(strcmp(domain.upperBoundary,'essential') || strcmp(domain.rightBoundary,'essential'))
    %if any of the two boundaries is essential, the node is essential
    nc(3,i) = 0;
    if(strcmp(domain.upperBoundary,'essential') && strcmp(domain.rightBoundary,'essential'))
    nc(4,i) = 0.5*(physical.Tb(3) + physical.Tb(2)); %We take the mean value of boundaries
    elseif(strcmp(domain.upperBoundary,'essential'))
    nc(4,i) = physical.Tb(3);
    elseif(strcmp(domain.rightBoundary,'essential'))
    nc(4,i) = physical.Tb(2);
    else
    error('No essential boundary?');
    end
    elseif(strcmp(domain.upperBoundary,'natural') && strcmp(domain.rightBoundary, 'natural'))
    %if both boundaries are natural, the node is natural
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = 0.5*(physical.qb(3) + physical.qb(2)); %we take the mean value
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(i == (domain.Nx + 1)*domain.Ny + 1) %upper left corner
    boundaryflag = true;
    if(strcmp(domain.upperBoundary,'essential') || strcmp(domain.leftBoundary,'essential'))
    %if any of the two boundaries is essential, the node is essential
    nc(3,i) = 0;
    if(strcmp(domain.upperBoundary,'essential') && strcmp(domain.leftBoundary,'essential'))
    nc(4,i) = 0.5*(physical.Tb(3) + physical.Tb(4)); %We take the mean value of boundaries
    elseif(strcmp(domain.upperBoundary,'essential'))
    nc(4,i) = physical.Tb(3);
    elseif(strcmp(domain.leftBoundary,'essential'))
    nc(4,i) = physical.Tb(4);
    else
    error('No essential boundary?');
    end
    elseif(strcmp(domain.upperBoundary,'natural') && strcmp(domain.leftBoundary, 'natural'))
    %if both boundaries are natural, the node is natural
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = 0.5*(physical.qb(3) + physical.qb(4)); %we take the mean value
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(i > 1 && i < domain.Nx + 1) %lower boundary
    boundaryflag = true;
    if(strcmp(domain.lowerBoundary,'essential'))
    nc(3,i) = 0;
    nc(4,i) = physical.Tb(1);
    elseif(strcmp(domain.lowerBoundary,'natural') && nc(3,i) ~= 0) %if nc(3,i) is already assigned 0, it is an 
                                                    %essential boundary and stays an essential boundary!
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = physical.qb(1);
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(mod(i,domain.Nx + 1) == 0 && i ~= domain.Nx + 1 && i ~= (domain.Nx + 1)*(domain.Ny + 1)) %right boundary
    boundaryflag = true;
    if(strcmp(domain.rightBoundary,'essential')) 
    nc(3,i) = 0;
    nc(4,i) = physical.Tb(2);
    elseif(strcmp(domain.rightBoundary,'natural') && nc(3,i) ~= 0)
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = physical.qb(2);
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end
end

if(i > (domain.Nx + 1)*domain.Ny + 1 && i < (domain.Nx + 1)*(domain.Ny + 1)) %upper boundary
    boundaryflag = true;
    if(strcmp(domain.upperBoundary,'essential'))
    nc(3,i) = 0;
    nc(4,i) = physical.Tb(3);
    elseif(strcmp(domain.upperBoundary,'natural') && nc(3,i) ~= 0)
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = physical.qb(3);
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(mod(i, domain.Nx + 1) == 1 && i ~= 1 && i ~= (domain.Nx + 1)*domain.Ny + 1) %left boundary
    boundaryflag = true;
    if(strcmp(domain.leftBoundary,'essential'))
    nc(3,i) = 0;
    nc(4,i) = physical.Tb(4);
    elseif(strcmp(domain.leftBoundary,'natural') && nc(3,i) ~= 0)
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
    nc(5,i) = physical.qb(4);
    elseif(isnan(nc(3,i)))
    error('Neither essential nor natural boundary');
    end 
end

if(~boundaryflag) %internal node
    %Assign equation number j
    nc(3,i) = j;
    j = j + 1;
end

    x = x + domain.elementLength;
    if(x - domain.length > p)
        x = 0;
        y = y + domain.elementLength;
    end
end

%single conversion for economic use of memory
nc = single(nc);

end

