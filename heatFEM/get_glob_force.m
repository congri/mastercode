function [F] = get_glob_force(domain, physical, k)
%Assemble global force vector

neq = max(domain.nodalCoordinates(3,:));
F = zeros(neq,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nested function for performance
Tb = zeros(4,1);
Tbflag = false;
function [f] = get_loc_force2(e, domain, kin, physical)
%Gives local force vector

%Contribution due to essential boundaries
%Boundary value temperature of element e
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
    Tb = zeros(4,1);
    Tbflag = false;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





for e = 1:domain.nElements
    f = get_loc_force2(e, domain, k, physical);
    if(sum(f))
        for ln = 1:4
           if(f(ln))
               eqn = domain.lm(e,ln);
               if(eqn)
                   F(eqn) = F(eqn) + f(ln);
               end
           end
        end
    end
end

end

