
a = rand(1, conductivity.dim);
[lambda] = computeLambda(a, domain, conductivity.lambdaCutoff);
    

for i = 1:200
    [logU, gradLogU, U] = cont_ref_output(lambda, a, conductivity, physical, domain);
end
