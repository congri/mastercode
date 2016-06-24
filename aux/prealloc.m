%preallocation script

muArray = zeros(optim(1).maxIterations,conductivity.dim);
muArray(nIterations,:) = conductivity.mu_a;
logUMean = 0;
q_mean = 0;
grad_q_mean = zeros(1,conductivity.dim);
gradLogUMean = zeros(1,conductivity.dim);
logUmeanArray = zeros(optim(1).maxIterations,1);
q_meanArray = zeros(optim(1).maxIterations,1);
grad_q_meanArray = zeros(optim(1).maxIterations, conductivity.dim);
gradLogUMeanArray = zeros(optim(1).maxIterations, conductivity.dim);