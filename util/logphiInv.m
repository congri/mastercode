function [u] = logphiInv(X, u0)
%Inverse of logphi function,
%logphi = log(normcdf(z)).

u = fzero(@(u)(logphi(u) - X), u0);


end

