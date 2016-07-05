
a = rand(1, conductivity.dim);

grad_kl = zeros(4, 4, domain.nElements);

tic;
for l = 1:conductivity.dim
    for k = 1:domain.nElements
        grad_kl(:,:,k) = exp(a(l) + domain.exponent(k,l)')*conductivity.localStiffnessOffset;
    end
end
t1 = toc


tic;
exponent = domain.exponent;
for l = 1:conductivity.dim
    for k = 1:domain.nElements
        grad_kl(:,:,k) = exp(a(l) + exponent(k,l)')*conductivity.localStiffnessOffset;
    end
end
t2 = toc

tic;
exponent = domain.exponent';
for l = 1:conductivity.dim
    for k = 1:domain.nElements
        grad_kl(:,:,k) = exp(a(l) + exponent(k,l))*conductivity.localStiffnessOffset;
    end
end
t3 = toc

tic;
exponent = domain.exponent';
off = conductivity.localStiffnessOffset;
for l = 1:conductivity.dim
    for k = 1:domain.nElements
        grad_kl(:,:,k) = exp(a(l) + exponent(k,l))*off;
    end
end
t4 = toc