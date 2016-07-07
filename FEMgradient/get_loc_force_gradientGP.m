function [f] = get_loc_force_gradientGP(domain, kGradient, e)
%Gives local force vector gradient, given that the gradient of the local
%stiffness matrix is plugged in
%This function is equal to get_loc_force up to performance improvements due
%to gradient computation

f = zeros(4,1);
Tb = zeros(4,1);
Tbflag = 0;
en = domain.essNode;
%Boundary value temperature of element e
for i = 1:4
    if(en(e, i))

        Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
        Tbflag = true;

    end
end


if(Tbflag)
    f = -kGradient*Tb;
end



%THIS MAY BE BUGGED!!!! USE ABOVE VERSION IF IT IS NOT WORKING PROPERLY
% f = zeros(4,domain.nElements);
% for k = 1:size(domain.essNodes,1)
%    
%     e = domain.essNodes(k, 1);
%     f(:,e) = - kGradient(:,:,e)*domain.Tb(:,e);
%     
% end

% if(norm(f - f2))
% f
% f2
% Tb
% domain.Tb
% pause
% end




end