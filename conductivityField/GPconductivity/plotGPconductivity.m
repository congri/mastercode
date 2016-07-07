%plot script for GP conductivities

clear all;
n_el = 15;
x = linspace(1/(2*n_el), 1 - 1/(2*n_el), n_el);
[X2, X1] = meshgrid(x);
x = [X1(:) X2(:)]';

log_sigma_f2 = 0;
log_l1 = -3.8606;
log_l2 = -4.1711;
l1 = sqrt(exp(log_l1))
l2 = sqrt(exp(log_l2))
params = [log_l1, log_l2, log_sigma_f2];

%squared distances in x and y direction
Xsq = sq_dist(x(1,:));
Ysq = sq_dist(x(2,:));

K = GPcov(Xsq, Ysq, 'ardSE', params);
condest(K)
m = GPmean(x, 'zero', 0);

figure
colormap hot
for i = 1:4
    
    [lambda, lambdaField] = GPconductivity(m, K);
    s(i) = subplot(2,2,i);
    set(s(i), 'fontsize', 18)
%     [~, ctr(i)] = contourf(X1, X2, lambdaField, 256, 'linecolor', 'none');
    ctr(i) = pcolor(X1, X2, lambdaField);
    set(ctr(i), 'EdgeColor', 'none');
    hold on
    plot([.1 .1 + l1], [.1 .1], 'w', 'linewidth', 2)
    plot([.1 .1], [.1 .1 + l2], 'w', 'linewidth', 2)
    hold off
    caxis([0, 15])


    xlabel('x')
    ylabel('y')
    zlabel('\lambda')
    c = colorbar;
    ylabel(c, '\lambda')
    axis square
    
end