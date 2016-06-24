function [Lambda] = plotLambda(mu, domain, conductivity, physical, l)
%plot heat conductivity field

f = figure('name','Heat conductivity field realizations');
for j = 1:4
    subplot(2,2,j)
    lambda(:,j) = computeLambda(mu(j,:), domain, conductivity.lambdaCutoff);
    for i = 1:(l^2)
       column = mod((i - 1),l) + 1;
       row = floor((i - 1)/l) + 1;

       Lambda(row, column) = log(lambda(i,j));
    end
    X = linspace(0,1,l);
    X = meshgrid(X);
    Y = X';
%     contourf(X,Y,Lambda,linspace(0,100,256),'LineColor','none');
    contourf(X,Y,Lambda,256,'LineColor','none');
%     caxis([-7, 4])
    set(f, 'Position', [10, 280, 1000, 840]);
%     title('heat conductivity field');
%     xlabel('x');
%     ylabel('y');
    cb = colorbar;
%     ylabel(cb,'log(\lambda)');
    set(gca,'FontSize',28) 
    colormap('cool')
    
    hold on
%     Sx = size(physical.x);
%     for i = 1:Sx(1)
%     pt = plot(physical.x(i,1),physical.x(i,2),'x');
%     set(pt,'LineWidth',3,'markersize',14,'color','white');
%     end
% xSqr = [physical.sqLo, physical.sqHi, physical.sqHi, physical.sqLo, physical.sqLo];
% ySqr = [physical.sqLo, physical.sqLo, physical.sqHi, physical.sqHi, physical.sqLo];
% ptl = plot(xSqr,ySqr);
% set(ptl,'LineWidth',3,'color','white');
    hold off
end

f = figure('name','Temperature field realizations');
for j = 1:4
    X = linspace(0,1,l + 1);
    X = meshgrid(X);
    Y = X';
    
    D = zeros(2,2,domain.nElements);
    for i = 1:domain.nElements
        D(:,:,i) = lambda(i,j)*eye(2); %only isotropic material
    end

    %Now the heat conductivities are properly set up
    control.plt = 0;
    control.gradientComputation = 1;
    FEMout = heat2d(domain, physical, control, D);

    subplot(2,2,j)
%     contourf(X,Y,FEMout.Tff,256,'LineColor','none')
    surf(X,Y,FEMout.Tff)
%     contourf(X,Y,heatMap.Tff,linspace(-600,1000,256),'LineColor','none')
%     caxis([0, 20])
    set(f, 'Position', [10, 280, 1000, 840]);
%     title('Temperature field');
%     xlabel('x');
%     ylabel('y');
    cb = colorbar;
%     ylabel(cb,'T');
    set(gca,'FontSize',28) 
    colormap('jet')
    
    hold on

    hold off
end


end

