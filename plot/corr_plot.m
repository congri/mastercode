%subplot script for correlated convergence plots
close all;
f = figure;
set(f, 'Position', [500, 300, 1200, 768]); 

subplot(2,2,1)
semilogy(covarNormArray(1:nIterations))
title('Convergence of (\Sigma_i:\Sigma_i)^{1/2}')
xlabel('iteration i')
ylabel('(\Sigma_i:\Sigma_i)^{1/2}')
xlim([0 nIterations])
% hold on
% semilogy(drop)
% hold off


subplot(2,2,2)
plot(muArray(1:nIterations,:))
title('Convergence of \mu_i')
xlabel('iteration i')
ylabel('\mu_i')
xlim([0 nIterations])


subplot(2,2,3)
plot(log(qMeanArray(1:nIterations)))
title('log of expected value of V(\theta_i) under q_i')
xlabel('iteration i')
ylabel('log(<V(\theta_i)>_{q_i})')
xlim([0 nIterations])


subplot(2,2,4)
errorbar(UMeanArray(1:nIterations), sqrt(varUArray(1:(nIterations))))
title('Expected value of output U under q_i')
xlabel('iteration i')
ylabel('<U>_{q_i}')
xlim([0 nIterations])
