%principal parameter file

rng('shuffle')

%Do not change order of params scripts! For a new FEM mesh each script
%after boundaryConditions has to be reexecuted

%optimization parameters go here
optimParams;
%specify boundary conditions
boundaryConditions;
%FEM domain specification
domain.nElements = 100; %number of finite elements
domainParams;
%specify target for optimization
target;
%some physical params (heat source field etc.)
physicalParams;
%params of the heat conductivity field
conductivityParams;
%second domain script has to be executed after conductivity field params
domainParams2;



