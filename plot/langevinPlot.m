%subplot script
close all;
f = figure;
set(f, 'Position', [500, 300, 1200, 768]); 

subplot(2,2,1)
semilogy(muDiffArray(1:nIterations))
title('|\mu_i - \mu_{i - 1}| ')
xlabel('iteration i')
ylabel('|\mu_i - \mu_{i - 1}|')
xlim([0 nIterations])

subplot(2,2,2)
plot(TMeanArray(:,1:nIterations)')
title('Temperature Measurements')
xlabel('iteration i')
ylabel('T')
xlim([0 nIterations])
hold on
plot([0 nIterations], [21 21], 'k', [0 nIterations], [15 15], 'k')
hold off



subplot(2,2,3)
plot(gradientMeanArray(1:(nIterations),:))
title('Gradient \nabla_a log q(a)')
xlabel('iteration i')
ylabel('<\nabla_a log q(a)>_{q_i}')
xlim([0 nIterations])


subplot(2,2,4)
plot(dlogU_daMeanArray(1:(nIterations),:))
title('\nabla_a log U(a)')
xlabel('iteration i')
ylabel('<\nabla_a log U(a)>_{q_i}')
xlim([0 nIterations])
