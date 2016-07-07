function [coord_mat] = get_el_coord(domain)
%Get coordinate of the center of an element given the global element number
%ATTENTION! ONLY WORKS FOR SQUARES!

%element vector
n_el = 1:domain.nElements;
%element number in x direction
nx = mod((n_el - 1),domain.Nx);
%element number in y direction
ny = floor((n_el - 1)/domain.Nx);
%coordinates
x = (nx + .5)*domain.elementLength;
y = (ny + .5)*domain.elementLength;

coord_mat = [x; y; n_el];

%convert to single for efficient memory usage
% coord_mat = single(coord_mat);


end