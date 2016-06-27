function [closestNodes] = getClosestNode(domain, physical)
%Computes the closest FEM nodes to the measurement points specified in
%physical.x

closestNodes = zeros(1,length(physical.T_target),'int16');
for i = 1:length(physical.T_target)
    xDiff = domain.nodalCoordinates(1,:) - physical.x(i,1);
    yDiff = domain.nodalCoordinates(2,:) - physical.x(i,2);
    nodalDistancesSquared = xDiff.^2 + yDiff.^2;
    [~,closestNodes(i)] = min(nodalDistancesSquared);   
end


end

