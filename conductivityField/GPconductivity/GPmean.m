function [m] = GPmean(x, type, params)
%Gives the Gaussian process mean function at x
    

n_el = size(x, 2);
if strcmp('zero', type)
    
    m = zeros(n_el, 1);
end
    
end

