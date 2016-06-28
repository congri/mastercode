function [f] = get_loc_force(e, domain, kin, physical)
%Gives local force vector

%Contribution due to essential boundaries

%Boundary value temperature of element e
Tb = zeros(4,1);
Tbflag = false;
for i = 1:4
    if(~isnan(domain.nodalCoordinates(4,domain.globalNodeNumber(e,i))))
        Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
        Tbflag = true;
    end
end

if(physical.source)
    f = physical.fh(:,e) + physical.fs(:,e);
else
    f = physical.fh(:,e);
end

if(Tbflag)
    fT = kin(:,:,e)*Tb;
    f = f - fT;
end

end

