%plot script for GP conductivities

clear all;
n_el = 10;
x = linspace(1/(2*n_el), 1 - 1/(2*n_el), n_el);
[X2, X1] = meshgrid(x);
x = [X1(:) X2(:)]';

log_sigma_f2 = 0;
log_l1 = -3.4004;
log_l2 = -3.6707;
params = [log_l1, log_l2, log_sigma_f2];

%squared distances in x and y direction
Xsq = sq_dist(x(1,:));
Ysq = sq_dist(x(2,:));

K = GPcov(Xsq, Ysq, 'ardSE', params);
rcond(K)
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

    xlabel('x')
    ylabel('y')
    zlabel('\lambda')
    c = colorbar;
    ylabel(c, '\lambda')
    axis square
    
end