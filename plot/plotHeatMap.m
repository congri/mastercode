%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot command
f = figure('name','2D heat map');
set(f, 'Position', [10, 350, 720, 540]);
contourf(X,Y,Tff,linspace(-50,200,256),'LineColor','none');
title('2D heat map');
xlabel('x');
ylabel('y');
cb = colorbar;
ylabel(cb,'T');
set(gca,'FontSize',14) 
colormap('jet')
dim = [.15 .6 .3 .3];
%concatenate textbox string
if(strcmp(domain.lowerBoundary,'essential'))
    s = 'T_{lo} =';
    s = [s ' ' num2str(physical.Tb(1))];
else
    s = 'q_{lo} = ';
    s = [s ' ' num2str(physical.qb(1))];
end
if(strcmp(domain.rightBoundary,'essential'))
    s = [s ' ' '\nT_{r} ='];
    s = [s ' ' num2str(physical.Tb(2))];
else
    s = [s ' ' '\nq_{r} = '];
    s = [s ' ' num2str(physical.qb(2))];
end
if(strcmp(domain.upperBoundary,'essential'))
    s = [s ' ' '\nT_{u} ='];
    s = [s ' ' num2str(physical.Tb(3))];
else
    s = [s ' ' '\nq_{u} = '];
    s = [s ' ' num2str(physical.qb(3))];
end
if(strcmp(domain.leftBoundary,'essential'))
    s = [s ' ' '\nT_{le} ='];
    s = [s ' ' num2str(physical.Tb(4))];
else
    s = [s ' ' '\nq_{le} = '];
    s = [s ' ' num2str(physical.qb(4))];
end
str = sprintf(s);
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold on
Sx = size(physical.x);
for i = 1:Sx(1)
pt = plot(physical.x(i,1),physical.x(i,2),'x');
set(pt,'LineWidth',3,'markersize',14,'color','white');
end
hold off