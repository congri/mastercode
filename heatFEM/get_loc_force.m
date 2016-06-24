function [f] = get_loc_force(e, domain, kin, physical)
%Gives local force vector

    %Contribution due to essential boundaries
    %local stiffness matrix k
    k = kin(:,:,e);
    
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
        fT = k*Tb;
        f = physical.fh(:,e) + physical.fs(:,e) - fT;
    else
        f = physical.fh(:,e) + physical.fs(:,e);
    end

end

