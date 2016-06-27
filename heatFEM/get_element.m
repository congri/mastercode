function [el] = get_element(domain, physical)
%plug in list of coordinates (x in first column, y in second column), get
%out corresponding list of element numbers of where x is located
%ONLY WORKS FOR SQUARE DOMAINS WITH SQUARE ELEMENTS!!!

%number of points in x
N_points = length(physical.T_target);

%Go through all points and determine elements
el = zeros(N_points,1,'int16');
for i = 1:N_points
   row = floor(physical.x(i,2)/domain.elementLength);
   column = ceil(physical.x(i,1)/domain.elementLength);
   
   el(i) = column + domain.Nx*row;  
end


end

