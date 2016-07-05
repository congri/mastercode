%plot script for GP conductivities

clear all;
n_el = 50;
x = linspace(1/(2*n_el), 1 - 1/(2*n_el), n_el);
[X2, X1] = meshgrid(x);
x = [X1(:) X2(:)]';

log_sigma_f2 = 0;
log_l1 = -5;
log_l2 = 0;
params = [log_l1, log_l2, log_sigma_f2];
k = @(x) GPcov(x, 'ardSE', params);
mu = @(x) GPmean(x, 'zero', 0);


figure
colormap hot
for i = 1:4
    
    [lambda, lambdaField] = GPconductivity(mu, k, x);
    s(i) = subplot(2,2,i);
    set(s(i), 'fontsize', 18)
    [~, ctr(i)] = contourf(X1, X2, lambdaField, 256, 'linecolor', 'none');
    xlabel('x')
    ylabel('y')
    zlabel('\lambda')
    c = colorbar;
    ylabel(c, '\lambda')
    axis square
    
end