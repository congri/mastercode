%plot shape function interpolate

% a = 4*rand(1, 9) - 3;
%compute conductivities
[lambda] = computeLambda(a, domain, conductivity.lambdaCutoff);
D = zeros(2,2,domain.nElements);
for i = 1:domain.nElements
    D(:,:,i) = lambda(i)*eye(2); %only isotropic material
end

%FEM
control.plt = 0;
control.gradientComputation = 0;
FEMout = heat2d(domain, physical, control, D);

%interpolation mesh
nint = 50;
d = 1e-6;
x = linspace(d, 1 - d, nint);
[X, Y] = meshgrid(x);

%element row and column of each coordinate
col = floor(X./domain.elementLength) + 1;
row = floor(Y./domain.elementLength) + 1;

el = col + domain.Nx*(row - 1);

%shape functions
A = domain.elementLength^2;
N1 = @(x, y, el) (1/A)*(x - domain.lc(el, 2, 1))*(y - domain.lc(el, 4, 2));
N2 = @(x, y, el) -(1/A)*(x - domain.lc(el, 1, 1))*(y - domain.lc(el, 4, 2));
N3 = @(x, y, el) (1/A)*(x - domain.lc(el, 1, 1))*(y - domain.lc(el, 1, 2));
N4 = @(x, y, el) -(1/A)*(x - domain.lc(el, 2, 1))*(y - domain.lc(el, 1, 2));


%shape function derivatives
B1x = @(x, y, el) (1/A)*(y - domain.lc(el, 4, 2));
B1y = @(x, y, el) (1/A)*(x - domain.lc(el, 2, 1));

B2x = @(x, y, el) -(1/A)*(y - domain.lc(el, 4, 2));
B2y = @(x, y, el) -(1/A)*(x - domain.lc(el, 1, 1));

B3x = @(x, y, el) (1/A)*(y - domain.lc(el, 1, 2));
B3y = @(x, y, el) (1/A)*(x - domain.lc(el, 1, 1));

B4x = @(x, y, el) -(1/A)*(y - domain.lc(el, 1, 2));
B4y = @(x, y, el) -(1/A)*(x - domain.lc(el, 2, 1));


%global node number array assembly
gnntemp = linspace(1, domain.Nx + 1, domain.Nx + 1);
for i = 1:(domain.Ny + 1)
    
    gnn(i, :) = gnntemp + (domain.Ny + 1 - i)*gnntemp(end);
    
end

[~, index] = sort(gnn(:));
%Tf takes global node number and gives temperature
Tf = FEMout.Tff(index);

Tfield = zeros(nint);
gradFieldX = zeros(nint);
gradFieldY = zeros(nint);
fluxX = zeros(nint);
fluxY = zeros(nint);
for i = 1:numel(X)

    T1 = Tf(domain.globalNodeNumber(el(i), 1));
    T2 = Tf(domain.globalNodeNumber(el(i), 2));
    T3 = Tf(domain.globalNodeNumber(el(i), 3));
    T4 = Tf(domain.globalNodeNumber(el(i), 4));

    Tfield(i) = Tfield(i) + T1*N1(X(i), Y(i), el(i)) + T2*N2(X(i), Y(i), el(i))...
                    + T3*N3(X(i), Y(i), el(i)) + T4*N4(X(i), Y(i), el(i));

    gradFieldX(i) = gradFieldX(i) + T1*B1x(X(i), Y(i), el(i)) + T2*B2x(X(i), Y(i), el(i))...
                    + T3*B3x(X(i), Y(i), el(i)) + T4*B4x(X(i), Y(i), el(i));
    gradFieldY(i) = gradFieldY(i) + T1*B1y(X(i), Y(i), el(i)) + T2*B2y(X(i), Y(i), el(i))...
                    + T3*B3y(X(i), Y(i), el(i)) + T4*B4y(X(i), Y(i), el(i));

end

for i = 1:numel(X)
    fluxX(i) = - lambda(el(i))*gradFieldX(i);
    fluxY(i) = - lambda(el(i))*gradFieldY(i);
end

figure
colormap jet
subplot(1,2,1)
contourf(X, Y, Tfield, 256, 'linecolor', 'none')
axis square
Xr = linspace(0, 1, domain.Nx + 1);
[Xr, Yr] = meshgrid(Xr);
s = subplot(1,2,2)
contourf(Xr, Yr, FEMout.Tff, 256, 'linecolor', 'none')
axis square
figure
subplot(1,2,1)
contour(X, Y, Tfield, 12)
hold on
quiver(X, Y, gradFieldX, gradFieldY, 'linewidth',1, 'markersize', 12)
hold off
xlim([0 1])
ylim([0 1])
axis square
subplot(1,2,2)
contour(X, Y, Tfield, 12)
hold on
quiver(X, Y, fluxX, fluxY, 'linewidth',1, 'markersize', 12)
hold off
xlim([0 1])
ylim([0 1])
axis square





