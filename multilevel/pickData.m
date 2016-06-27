function [approx_uOut, inputDataOut] = pickData(approx_uIn, inputDataIn, nData)
%Pick out equally spaced data

hi = max(approx_uIn);
lo = min(approx_uIn);

p = .9;
centers = linspace(p*hi + (1 - p)*lo,hi,nData);

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

