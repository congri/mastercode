function [k, grad] = GPcov(X2, Y2, type, params)
%Gives GP covariance for heat conductivity field
%   X:  Mutual squared distances in x-direction of finite element centers
%   Y:  Mutual squared distances in y-direction of finite element centers

if strcmp(type, 'isoSE')
    
    log_sigma_f2 = params(2);
    l2 = exp(params(1));
    
    k = exp(log_sigma_f2 - X2/l2 - Y2/l2);

    if(nargout > 1) %compute derivatives

        %w.r.t. l^2
        dk(:,:,1) = (X2/(l2^2) + Y2/(l2^2)).*k;
        %w.r.t. sigma_f^2
        dk(:,:,2) = exp(- X2/l2 - Y2/l2);

    end
    
elseif strcmp(type, 'ardSE')
    
    log_sigma_f2 = params(3);
    l2 = exp(params(1:2));
    
    k = exp(log_sigma_f2 - X2/l2(1) - Y2/l2(2));

    log_threshold = -25;
    k(log_sigma_f2 - X2/l2(1) - Y2/l2(2) < log_threshold) = 0;
    k = sparse(k);

    if(nargout > 1) %compute derivatives

%         %w.r.t. l_1^2
%         dk(:,:,1) = X2/(l2(1)^2).*k;
%         %w.r.t. l_2^2
%         dk(:,:,2) = Y2/(l2(2)^2).*k;
%         %w.r.t. sigma_f^2
%         dk(:,:,3) = exp(- X2/l2(1) - Y2/l2(2));


        %w.r.t. l_1^2
        grad(1).dK = X2/(l2(1)^2).*k;
        %w.r.t. l_2^2
        grad(2).dK = Y2/(l2(2)^2).*k;
        %w.r.t. sigma_f^2
        grad(3).dK = exp(- X2/l2(1) - Y2/l2(2));

    end
    
end
    
    
end

