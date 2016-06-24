function [logU, transformedU] = generateData(inputData, physical, domain, conductivity, nData, namestr)
%Generates data samples for multilevel code

%prealloc
logU = zeros(nData,1);
transformedU = zeros(nData,1);

for i = 1:nData
    i
    lambda = computeLambda(inputData(i,:), domain, conductivity.lambdaCutoff);

    %compute output
    [logUi] = cont_ref_output(lambda, inputData(i,:), conductivity, physical, domain);
    logU(i) = logUi;
    transformedU(i) = logphiInv(logUi,-20);
end

%save data
mkdir('./data/multilevel')
namestr = strcat('./data/multilevel/',namestr);
save(namestr, 'logU', 'transformedU','-v7.3');

end

