function [Phi] = designMatrix(dataIn, nData)
%Compute the design matrix with linear basis functions

Phi = ones(nData, 2);
Phi(:,2) = dataIn;


end

