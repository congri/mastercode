function [approx_uOut, inputDataOut] = pickData(approx_uIn, inputDataIn, nData, upperFraction)
%Pick out equally spaced data

[approx_uIn, index] = sort(approx_uIn);
inputDataIn = inputDataIn(index,:);
% lu = length(approx_uIn);
% approx_uIn = approx_uIn((1 - upperFraction)*lu:end);
% inputDataIn = inputDataIn((1 - upperFraction)*lu:end,:);

hi = max(approx_uIn);
lo = min(approx_uIn);

% p = 1 - upperFraction;
centers = linspace(lo,hi,nData);

approx_uOut = zeros(nData,1);
inputDataOut = zeros(nData, size(inputDataIn,2));
for i = 1:nData
    
   [~,index] = min(abs(centers(i) - approx_uIn));
   approx_uOut(i) = approx_uIn(index);
   approx_uIn(index) = [];      %each data point can only be taken once
   inputDataOut(i,:) = inputDataIn(index,:);
   inputDataIn(index,:) = [];
   
end

end

