
a = rand(1, conductivity.dim);

grad_kl = zeros(4, 4, domain.nElements);

% tic;
% for l = 1:conductivity.dim
%     for k = 1:domain.nElements
%         grad_kl(:,:,k) = exp(a(l) + domain.exponent(k,l)')*conductivity.localStiffnessOffset;
%     end
% end
% t1 = toc
% 
% 
% tic;
% exponent = domain.exponent;
% for l = 1:conductivity.dim
%     for k = 1:domain.nElements
%         grad_kl(:,:,k) = exp(a(l) + exponent(k,l)')*conductivity.localStiffnessOffset;
%     end
% end
% t2 = toc
% 
% tic;
% exponent = domain.exponent';
% for l = 1:conductivity.dim
%     for k = 1:domain.nElements
%         grad_kl(:,:,k) = exp(a(l) + exponent(k,l))*conductivity.localStiffnessOffset;
%     end
% end
% t3 = toc
% 
% tic;
% exponent = domain.exponent';
% off = conductivity.localStiffnessOffset;
% for l = 1:conductivity.dim
%     for k = 1:domain.nElements
%         grad_kl(:,:,k) = exp(a(l) + exponent(k,l))*off;
%     end
% end
% t4 = toc


tic;
for j = 1:1
    f = zeros(4,domain.nElements);
    Tb = zeros(4,1);
    Tbflag = 0;
    for e = 1:domain.nElements
        %Boundary value temperature of element e
        for i = 1:4
            %             if(~isnan(domain.nodalCoordinates(4,domain.globalNodeNumber(e,i))))
            if(domain.essNode(e,i))
                Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
                Tbflag = true;
            end
        end
        
        
        if(Tbflag)
            f(:,e) = -grad_kl(:,:,e)*Tb;
            Tb = zeros(4,1);
            Tbflag = false;
        end
    end
end
t5 = toc


tic;
for j = 1:1
    f = zeros(4,domain.nElements);
    Tb = zeros(4,1);
    ess = domain.essNode;
    gN = domain.globalNodeNumber;
    Tbflag = 0;
    for e = 1:domain.nElements
        %Boundary value temperature of element e
        for i = 1:4
            %             if(~isnan(domain.nodalCoordinates(4,domain.globalNodeNumber(e,i))))
            if(ess(e, i))
                Tb(i) = domain.nodalCoordinates(4,gN(e,i));
                Tbflag = true;
            end
        end
        
        
        if(Tbflag)
            f(:,e) = -grad_kl(:,:,e)*Tb;
            Tb = zeros(4,1);
            Tbflag = false;
        end
    end
end
t6 = toc


k = 1;
for e = 1:domain.nElements
    %Boundary value temperature of element e
    for i = 1:4
        
        if(domain.essNode(e,i))
%             Tb(i) = domain.nodalCoordinates(4,domain.globalNodeNumber(e,i));
%             Tbflag = true;

            E(k) = e;
            I(k) = i;
            
            k = k + 1;
        end
        
    end
    
    
    if(Tbflag)
        f(:,e) = -grad_kl(:,:,e)*Tb;
        Tb = zeros(4,1);
        Tbflag = false;
    end
end