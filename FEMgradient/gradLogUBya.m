function [dlogU_da] = gradLogUBya(FEMout, conductivity, domain, physical, a)
%Computes the gradient d log U(a)/d a
[gradK, gradF] = gradKgradF(FEMout, conductivity, domain, a);


%compute (d/dT) log U(T), where T are natural node temperatures
NEquations = length(FEMout.naturalTemperatures);
%derivative of measurement temperature by natural node temperatures
dTMeasure_dTNatNode = zeros(length(physical.T_target), NEquations);
for i = 1:physical.numberOfMeasurements
    if(domain.nodalCoordinates(3,physical.closestNodes(i)) > 0)
        dTMeasure_dTNatNode(i, domain.nodalCoordinates(3,physical.closestNodes(i))) = 1;
    end
end

dLogU_dT = -(FEMout.temperatureMeasurements - physical.T_target)'*physical.covTargetInv*dTMeasure_dTNatNode;

%Compute adjoints
adjoints = FEMout.globalStiffness\dLogU_dT'; %adjoints seem to be correct 18/02/16

%compute d log U(a)/d a
dlogU_da = zeros(1,conductivity.dim);
for l = 1:conductivity.dim
    dlogU_da(l) = -adjoints'*(gradK(:,:,l)*FEMout.naturalTemperatures - gradF(:,l));
end    



end

