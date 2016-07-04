function [gfd, g, rel] = gradcheck(f, x, delta)
%finite difference gradient check. function handle f is a handle to a
%function giving [f, df], i.e. function value and gradient


[fun, g] = f(x);

lx = length(x);
dx = zeros(1, lx);
gfd = zeros(1, lx);
rel = zeros(1, lx);

for i = 1:lx
   
    dx(i) = delta;
    fun2 = f(x + dx);
    gfd(i) = (fun2 - fun)/delta;
    dx(i) = 0;
    
end
    

rel = gfd./g;
end

