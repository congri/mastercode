

tic
for a = 1:500
f = zeros(4,domain.nElements);
Tb = zeros(4,1);
Tbflag = 0;
for e = 1:domain.nElements
    %Boundary value temperature of element e
    for i = 1:4
        if(domain.essNode(e,i))
            
            Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
            Tbflag = true;
            
        end
    end
    
    
    if(Tbflag)
        f(:,e) = -kGradient(:,:,e)*Tb;
        Tb = zeros(4,1);
        Tbflag = false;
    end
end
end
t1 = toc


tic
for a = 1:500
f = zeros(4,domain.nElements);
Tb = zeros(4,1);
Tbflag = 0;
en = domain.essNode;
for e = 1:domain.nElements
    %Boundary value temperature of element e
    for i = 1:4
        if(en(e, i))
            
            Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
            Tbflag = true;
            
        end
    end
    
    
    if(Tbflag)
        f(:,e) = -kGradient(:,:,e)*Tb;
        Tb = zeros(4,1);
        Tbflag = false;
    end
end

t2 = toc