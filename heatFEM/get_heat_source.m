function [fs] = get_heat_source(s,domain)
%Gets the elements of the local force due to the heat source (an array with
%input element number e and local node number i

%Edge length of element
Ae = domain.elementLength^2;

%Gauss points
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
eta1 = -1/sqrt(3);
eta2 = 1/sqrt(3);

fs = zeros(4,domain.nElements,'single');

for e = 1:domain.nElements
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


    fs(1,e) = s(e)*(1/Ae)*((xI - x2)*(yI - y4) + (xII - x2)*(yII - y4) + (xI - x2)*(yII - y4) + (xII - x2)*(yI - y4));
    fs(2,e) = -s(e)*(1/Ae)*((xI - x1)*(yI - y4) + (xII - x1)*(yII - y4) + (xI - x1)*(yII - y4) + (xII - x1)*(yI - y4));
    fs(3,e) = s(e)*(1/Ae)*((xI - x1)*(yI - y1) + (xII - x1)*(yII - y1) + (xI - x1)*(yII - y1) + (xII - x1)*(yI - y1));
    fs(4,e) = -s(e)*(1/Ae)*((xI - x2)*(yI - y1) + (xII - x2)*(yII - y1) + (xI - x2)*(yII - y1) + (xII - x2)*(yI - y1));
end

end

