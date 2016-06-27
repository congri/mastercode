function [Equations, LocalNode] = get_equations(nElements, lm)
%Equation number array for sparse global stiffness assembly


%Can this be done more efficiently?
%Assign global stiffness matrix
localNodeInit = 1:4;
Equations = [];
LocalNode = [];
for e = 1:nElements
    equations = lm(e,localNodeInit);
    equations = equations(equations > 0);
    localNode = localNodeInit(equations > 0);

    [Equations1, Equations2] = meshgrid(equations);
    Equations = [Equations; Equations1(:) Equations2(:)];

    [LocalNode1, LocalNode2] = meshgrid(localNode);
    LocalNode = [LocalNode; LocalNode1(:) LocalNode2(:) repmat(e, length(equations)^2, 1)];
end

end

