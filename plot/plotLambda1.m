function [Lambda] = plotLambda1(mu, domain, conductivity, physical, l)
%plot heat conductivity field

f = figure('name','Heat conductivity field realizations');
% sp1 = subplot(2,1,1)
lambda = computeLambda(mu, domain, conductivity.lambdaCutoff);
for i = 1:(l^2)
   column = mod((i - 1),l) + 1;
   row = floor((i - 1)/l) + 1;

   Lambda(row, column) = log(lambda(i));
end
X = linspace(0,1,l);
X = meshgrid(X);
Y = X';
%     contourf(X,Y,Lambda,linspace(0,100,256),'LineColor','none');
[~,ctf1] = contourf(X,Y,Lambda,256,'LineColor','none');
%     caxis([-7, 5])
set(f, 'Position', [10, 350, 720, 540]);
title('heat conductivity field');
xlabel('x');
ylabel('y');
cb1 = colorbar;
ylabel(cb1,'log(\lambda)');
set(gca,'FontSize',24) 
colormap('cool')

hold on
% Sx = size(physical.x);
% for i = 1:Sx(1)
% pt = plot(physical.x(i,1),physical.x(i,2),'x');
% plot(linspace(min(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),min(physical.x(:,2)),10),'r')
% set(pt,'LineWidth',3,'markersize',14,'color','white');
% end
%plot middle plaquette
% ptl = plot(linspace(min(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),min(physical.x(:,2)),10),...
%     linspace(min(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(max(physical.x(:,2)),max(physical.x(:,2)),10),...
%     linspace(min(physical.x(:,1)),min(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),max(physical.x(:,2)),10),...
%     linspace(max(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),max(physical.x(:,2)),10));
% set(ptl,'LineWidth',3,'color','white');
hold off


f = figure('name','Temperature field realizations');
% subplot(2,1,2)
X = linspace(0,1,l + 1);
X = meshgrid(X);
Y = X';
heatMap = cont_ref_output(lambda, zeros(1,size(mu,2)), conductivity, physical, domain);
max(max(heatMap.Tff))
min(min(heatMap.Tff))
% heatMap.temperatureMeasurements
heatMap.Tff(heatMap.Tff > 50) = 50;
heatMap.Tff(heatMap.Tff < 5) = 5;
contourf(X,Y,heatMap.Tff,256,'LineColor','none')
%     contourf(X,Y,heatMap.Tff,linspace(-600,1000,256),'LineColor','none')
%  caxis([5, 50])
set(f, 'Position', [10, 350, 720, 540]);
title('Temperature field');
xlabel('x');
ylabel('y');
cb = colorbar;
ylabel(cb,'T');
set(gca,'FontSize',24) 
colormap('jet')

hold on
% Sx = size(physical.x);
% for i = 1:Sx(1)
% pt = plot(physical.x(i,1),physical.x(i,2),'x');
% set(pt,'LineWidth',3,'markersize',14,'color','white');
% end
%plot middle plaquette
% ptl = plot(linspace(min(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),min(physical.x(:,2)),10),...
%     linspace(min(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(max(physical.x(:,2)),max(physical.x(:,2)),10),...
%     linspace(min(physical.x(:,1)),min(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),max(physical.x(:,2)),10),...
%     linspace(max(physical.x(:,1)),max(physical.x(:,1)),10)...
%     ,linspace(min(physical.x(:,2)),max(physical.x(:,2)),10));
% set(ptl,'LineWidth',3,'color','white');
hold off


end



