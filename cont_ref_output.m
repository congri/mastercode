function [logU, gradLogU, U] = cont_ref_output(lambda, aTemp, conductivity, physical, domain)
%Gives the output function in the continuous lambda case
%Construct heat conductivity matrix D

D = zeros(2,2,domain.nElements);
for i = 1:domain.nElements
    D(:,:,i) = lambda(i)*eye(2); %only isotropic material
end

%Now the heat conductivities are properly set up
control.plt = 0;
control.gradientComputation = 1;
FEMout = heat2d(domain, physical, control, D);
temperatureMeasurements = FEMout.temperatureMeasurements;

%Compute exponent explicitly
TDiff = (temperatureMeasurements - physical.T_target);
logU = -.5*TDiff'*physical.covTargetInv*TDiff;
if(nargout > 2)
    U = exp(logU);
end

%Gradient computation
if(nargout > 1)
    gradLogU = gradLogUBya(FEMout, conductivity, domain, physical, aTemp');
end

%Finite differences gradient check
gradcheck = 1;
if(gradcheck)
    disp('Perform finite difference gradient check')
    delta_a = 1e-4;
    FDgrad_dlogU_da = zeros(1,conductivity.dim);
    for d = 1:conductivity.dim
        deltaVec = zeros(1,conductivity.dim);
        deltaVec(d) = 1;
        lambdaTest = computeLambda(aTemp + delta_a*deltaVec, domain, conductivity.lambdaCutoff);
        DTest = zeros(2,2,domain.nElements);
        for i = 1:domain.nElements
            DTest(:,:,i) = lambdaTest(i)*eye(2); %only isotropic material
        end
        FEMoutTest = heat2d(domain, physical, control, DTest);
%         testOutput = physical.normalizationConstant*mvnpdf(FEMoutTest.temperatureMeasurements, physical.T_target, physical.cov_target);
        TDifftest = (FEMoutTest.temperatureMeasurements - physical.T_target);
        testLogU = -.5*TDifftest'*physical.covTargetInv*TDifftest;
        trueLogU = logU;
        outputDiff = testLogU - trueLogU;
        FDgrad_dlogU_da(d) = outputDiff/delta_a;
    end
    FDgrad_dlogU_da
    lambdaDiff = norm(lambda - lambdaTest)
    trueGrad = gradLogU
    relgrad = FDgrad_dlogU_da./trueGrad
    disp('----------------------------------------------------------')
end

end