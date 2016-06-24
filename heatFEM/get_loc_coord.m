function [lc] = get_loc_coord(d)
%Gives arrays taking the element and local node number and
%giving the nodal coordinate

lc = zeros(d.nElements,4,2);
for e = 1:d.nElements
        %x coordinates
        lc(e,1,1) = (e-1)*d.elementLength - floor((e - 1)/d.Nx)*d.Nx*d.elementLength;
        lc(e,2,1) = e*d.elementLength - floor((e - 1)/d.Nx)*d.Nx*d.elementLength;
        lc(e,3,1) = e*d.elementLength - floor((e - 1)/d.Nx)*d.Nx*d.elementLength;
        lc(e,4,1) = (e-1)*d.elementLength - floor((e - 1)/d.Nx)*d.Nx*d.elementLength;

        %y coordinates
        lc(e,1,2) = floor((e - 1)/d.Nx)*d.elementLength;
        lc(e,2,2) = floor((e - 1)/d.Nx)*d.elementLength;
        lc(e,3,2) = floor((e - 1)/d.Nx)*d.elementLength + d.elementLength;
        lc(e,4,2) = floor((e - 1)/d.Nx)*d.elementLength + d.elementLength;
end

end

