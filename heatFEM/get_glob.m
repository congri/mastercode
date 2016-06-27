function [g] = get_glob(d)
%Get global node number from global element number
%and local node number

g = zeros(d.nElements, 4);
for e = 1:d.nElements
    for l = 1:4
        g(e,1) = e + floor((e - 1)/d.Nx);
        g(e,2) = e + floor((e - 1)/d.Nx) + 1;
        g(e,3) = g(e,1) + d.Nx + 2;
        g(e,4) = g(e,1) + d.Nx + 1;
    end
end
    
end

