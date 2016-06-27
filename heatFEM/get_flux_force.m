function [fh] = get_flux_force(domain, physical)
%Contribution to local force due to heat flux

fh = sparse(4,domain.nElements);

%Gauss points
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
eta1 = -1/sqrt(3);
eta2 = 1/sqrt(3);

for e = 1:domain.nElements
    %Contribution due to free boundaries
    
    %short hand notation. Coordinates of local nodes
    x1 = domain.lc(e,1,1);
    x2 = domain.lc(e,2,1);
    y1 = domain.lc(e,1,2);
    y4 = domain.lc(e,4,2);
    
    %Coordinate transformation
    xI = 0.5*(x1 + x2) + 0.5*xi1*(x2 - x1);
    xII = 0.5*(x1 + x2) + 0.5*xi2*(x2 - x1);
    yI = 0.5*(y1 + y4) + 0.5*eta1*(y4 - y1);
    yII = 0.5*(y1 + y4) + 0.5*eta2*(y4 - y1);
    
    %Check if element belongs to natural boundary
    j = 0;
    for i = 1:4
        if(~isnan(domain.nodalCoordinates(5,domain.globalNodeNumber(e,i)))) %if it is a number, a heat flux is assigned --> natural boundary          
            j = 1;
            break;
        end
    end

    if(j > 0) %The element has natural boundaries

        if(e == 1) %lower left corner
           if(strcmp(domain.lowerBoundary,'natural'))
            h = physical.qb(1);
            fh(1,e) = fh(1,e) + h*(1/(2*l))*(-y4)*(xI + xII - 2*x2);
            fh(2,e) = fh(2,e) + h*(1/(2*l))*y4*(xI + xII - 2*x1);
           elseif(strcmp(domain.leftBoundary, 'natural'))
            h = physical.qb(4);
            fh(1,e) = fh(1,e) - h*(1/2)*(yI + yII - 2*y4);
            fh(4,e) = fh(4,e) + h*(1/2)*(yI + yII - 2*y1);
           else
               error('Which boundary of corner element is natural?');
           end

        elseif(e == domain.Nx) %lower right corner
            if(strcmp(domain.lowerBoundary,'natural'))
            h = physical.qb(1);
            fh(1,e) = fh(1,e) + h*(1/(2*l))*(-y4)*(xI + xII - 2*x2);
            fh(2,e) = fh(2,e) + h*(1/(2*l))*y4*(xI + xII - 2*x1);
           elseif(strcmp(domain.rightBoundary, 'natural'))
            h = physical.qb(2);
            fh(2,e) = fh(2,e) - h*(1/2)*(yI + yII - 2*y4);
            fh(3,e) = fh(3,e) + h*(1/2)*(yI + yII - 2*y1);
           else
               error('Which boundary of corner element is natural?');
            end

        elseif(e == domain.Nx*domain.Ny) %upper right corner
            if(strcmp(domain.upperBoundary,'natural'))
            h = physical.qb(3);
            fh(3,e) = fh(3,e) + h*(1/2)*(xI + xII - 2*x1);
            fh(4,e) = fh(4,e) - h*(1/2)*(xI + xII - 2*x2);
           elseif(strcmp(domain.rightBoundary, 'natural'))
            h = physical.qb(2);
            fh(2,e) = fh(2,e) - h*(1/2)*(yI + yII - 2*y4);
            fh(3,e) = fh(3,e) + h*(1/2)*(yI + yII - 2*y1);
           else
               error('Which boundary of corner element is natural?');
            end

        elseif(e == domain.Nx*(domain.Ny - 1) + 1) %upper left corner
            if(strcmp(domain.upperBoundary,'natural'))
            h = physical.qb(3);
            fh(3,e) = fh(3,e) + h*(1/2)*(xI + xII - 2*x1);
            fh(4,e) = fh(4,e) - h*(1/2)*(xI + xII - 2*x2);
           elseif(strcmp(domain.leftBoundary, 'natural'))
            h = physical.qb(4);
            fh(1,e) = fh(1,e) - h*(1/2)*(yI + yII - 2*y4);
            fh(4,e) = fh(4,e) + h*(1/2)*(yI + yII - 2*y1);
           else
               error('Which boundary of corner element is natural?');
            end

        elseif(e > 1 && e < domain.Nx && strcmp(domain.lowerBoundary,'natural')) %element is on lower boundary
            h = physical.qb(1);
            %Second order Gauss quadrature for linear function?
            fh(1,e) = fh(1,e) + h*(1/(2*l))*(-y4)*(xI + xII - 2*x2);
            fh(2,e) = fh(2,e) + h*(1/(2*l))*y4*(xI + xII - 2*x1);
        elseif(mod(e, domain.Nx) == 0 && e ~= domain.Nx && e ~= domain.Nx*domain.Ny && strcmp(domain.rightBoundary,'natural')) %element is on right boundary
            h = physical.qb(2);
            fh(2,e) = fh(2,e) - h*(1/2)*(yI + yII - 2*y4);
            fh(3,e) = fh(3,e) + h*(1/2)*(yI + yII - 2*y1);
        elseif(e > domain.Nx*(domain.Ny - 1) + 1 && e < domain.Nx*domain.Ny && strcmp(domain.upperBoundary,'natural')) %element is on upper boundary
            h = physical.qb(3);
            fh(3,e) = fh(3,e) + h*(1/2)*(xI + xII - 2*x1);
            fh(4,e) = fh(4,e) - h*(1/2)*(xI + xII - 2*x2);
        elseif(mod(e,domain.Nx) == 1 && e ~= 1 && e ~= domain.Nx*(domain.Ny - 1) + 1 && strcmp(domain.leftBoundary,'natural')) %element is on left boundary
            h = physical.qb(4);
            fh(1,e) = fh(1,e) - h*(1/2)*(yI + yII - 2*y4);
            fh(4,e) = fh(4,e) + h*(1/2)*(yI + yII - 2*y1);
        else
            error('Between which nodes is this boundary?'); 
        end
    end
end

end

