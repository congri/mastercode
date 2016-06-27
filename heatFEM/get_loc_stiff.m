function [k] = get_loc_stiff(e,lc,D)
%Gives the local stiffness matrix

%short hand notation
x1 = lc(e,1,1);
x2 = lc(e,2,1);
y1 = lc(e,1,2);
y4 = lc(e,4,2);

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
Z = zeros(2,2);
Dmat = [D Z Z Z; Z D Z Z; Z Z D Z; Z Z Z D];

k = Bvec'*Dmat*Bvec;
end

