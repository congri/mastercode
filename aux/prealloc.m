%preallocation script

muArray = zeros(optim(1).maxIterations,conductivity.dim);
muArray(nIterations,:) = conductivity.mu_a;
