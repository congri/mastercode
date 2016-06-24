function [Qx, Qy] = plotflux(X,Y,Tff,l,Nx,Ny,D)
%Plots the heat flux as quiver plot and gives the values back in matrices
%Qx, Qy


Qx = zeros(Nx,Ny);
Qy = zeros(Nx,Ny);
Qxt = zeros(Nx,Ny);
Qyt = zeros(Nx,Ny);
for i = 1:Nx
   for j = 1:Ny
      Qxt(i,j) = Tff(j,i+1) - Tff(j,i);
      Qyt(i,j) = Tff(j+1,i) - Tff(j,i);
      
      Qx(i,j) = -D(1,:,i + (j-1)*Nx)*[Qxt(i,j); Qyt(i,j)];
      Qy(i,j) = -D(2,:,i + (j-1)*Nx)*[Qxt(i,j); Qyt(i,j)]; 
   end
end

Qx = Qx'/l;
Qy = Qy'/l;

X(1,:) = [];
X(:,1) = [];
Y(1,:) = [];
Y(:,1) = [];

f = figure('name','Heat flux q');
set(f, 'Position', [735, 350, 720, 540]);
quiver(X,Y,Qx,Qy);
title('Heat flux q');
xlabel('x');
ylabel('y');
set(gca,'FontSize',14) 
xlim([0 1]);
ylim([0 1]);


end

