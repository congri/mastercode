function [f] = get_loc_force_gradient(domain, kGradient)
%Gives local force vector gradient, given that the gradient of the local
%stiffness matrix is plugged in
%This function is equal to get_loc_force up to performance improvements due
%to gradient computation

    %Contribution due to essential boundaries
    %local stiffness matrix k
    %k = kGradient(:,:,e);
    
    f = zeros(4,domain.nElements);
    for e = 1:domain.nElements
        %Boundary value temperature of element e
        Tb = zeros(4,1);
        Tbflag = 0;
        for i = 1:4
            if(~isnan(domain.nodalCoordinates(4,domain.globalNodeNumber(e,i))))
                Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
                Tbflag = 1;
            end
        end


        if(Tbflag)
            f(:,e) = -kGradient(:,:,e)*Tb;
        end
    end

end