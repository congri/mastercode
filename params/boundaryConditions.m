%specify boundary conditions here

physical.Tb = [0 0 0 0];
physical.qb = [.003 -.003 -.003 .001];

domain.lowerBoundary = 'essential';
domain.rightBoundary = 'natural';
domain.upperBoundary = 'essential';
domain.leftBoundary  = 'natural';