function [k] = GPcov(x, type, params)
%Gives GP covariance for heat conductivity field
%   x:  points where to draw the heat conductivity (i.e. centers of finite
%   elements

n_el = size(x, 2);
k = zeros(n_el);

if strcmp(type, 'isoSE')
    
    sigma_f2 = exp(params(2));
    l2 = exp(params(1));
    
    for i = 1:n_el
        for j = i:n_el
            
            k(i, j) = sigma_f2*exp(-(norm(x(:, i) - x(:, j))^2)/l2);
            k(j, i) = k(i, j);
            
        end
    end
    
elseif strcmp(type, 'ardSE')
    
    sigma_f2 = exp(params(3));
    l2 = exp(params(1:2));
    
    for i = 1:n_el
        for j = i:n_el
            
            k(i, j) = sigma_f2*exp(-((x(1, i) - x(1, j))^2)/l2(1) - ((x(2, i) - x(2, j))^2)/l2(2));
            k(j, i) = k(i, j);
            
        end
    end
    
end
    
    
end

