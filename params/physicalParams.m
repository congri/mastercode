%some physical params (heat source etc.)

%Assign heat source field
physical.heatSourceField = sparse(domain.nElements,1);
%Force contributions due to heat flux and source
physical.fs = get_heat_source(physical.heatSourceField, domain);
if(sum(sum(physical.fs)))
    physical.source = true;
else
    physical.source = false;
end
physical.fh = get_flux_force(domain, physical);
%Element numbers of measurement points
physical.measurementElements = get_element(domain, physical);
physical.closestNodes = getClosestNode(domain, physical);
physical.numberOfMeasurements = length(physical.T_target);