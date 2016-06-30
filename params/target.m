%target value specification to optimize on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %target middle square
% sqLo = .38;
% sqHi = .62;
% x = linspace(sqLo,sqHi,5);
% x = meshgrid(x);
% y = x';
% x = x(:);
% y = y(:);
% physical.x = [x, y; .05 .1; .05 .35; .05 .65; .05 .9; .95 .1; .95 .35; .95 .65; .95 .9];
% 
% clear sqLo sqHi x y;
% 
% physical.T_target = [10*ones(25,1); 15*ones(8,1)];
% physical.cov_target = (1/4000)*eye(length(physical.T_target));
% %regularize corners
% for i = 0:7
%     physical.cov_target((end - i),(end - i)) = 2500;
% end
% clear i;
% physical.covTargetInv = inv(physical.cov_target);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%target
[X,Y] = meshgrid(linspace(0,domain.length,domain.Nx + 1), linspace(0, domain.height, domain.Ny + 1));
physical.x = [X(:) Y(:)];
physical.x = single(physical.x);
% physical.T_target = 2*peaks(5*(X(:) - .5),5*(Y(:) - .5));
physical.T_target = 20*exp(-12*(5*X(:) + .5).*(Y(:) - .5).^2);
physical.cov_target = sparse((1/1)*eye(length(physical.T_target)));
physical.covTargetInv = inv(physical.cov_target);

