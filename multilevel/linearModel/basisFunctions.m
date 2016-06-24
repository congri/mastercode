function [f, df] = basisFunctions(input, options)
%Basis functions of the Bayesian linear model between multilevel code

f = zeros(options.nBasis, length(input));
df = zeros(options.nBasis, length(input));
if(strcmp(options.type,'polynomial'))
    %polynomial basis functions

    for i = 0:(options.nBasis - 1)
        f(i + 1,:) = input.^i;
        df(i + 1,:) = i*input.^(i - 1);
    end   
    
elseif(strcmp(options.type,'oddpolynomial'))
    %only odd polynomials as basis functions; ensures that u_s --> +- inf if u_f --> +- inf 
    
    for i = 0:(options.nBasis - 1)
        f(i + 1,:) = input.^(2*i + 1);
        df(i + 1,:) = (2*i + 1)*input.^(2*i);
    end  
    
elseif(strcmp(options.type,'hermite'))
    %hermite polynomial (physisist's)
    %!!!appears to be very slow!!! should be changed to explicit
    %polynomials
    
    for i = 0:(options.nBasis - 1)
        f(i + 1,:) = hermiteH(i,input);
        if(i > 0)
            df(i + 1,:) = 2*i*hermiteH(i - 1,input);
        else
            df(i + 1,:) = 0;
        end
    end 
    
elseif(strcmp(options.type,'rbf'))
    %Gaussian radial basis functions
    
    for i = 0:(options.nBasis - 1)
        f(i + 1,:) = exp(-(1/(2*options.rbfStd^2))*(input - options.rbfLocations(i + 1)).^2);
        df(i + 1,:) = -(1/options.rbfStd^2)*(input - options.rbfLocations(i + 1))*...
            exp(-(1/(2*options.rbfStd^2))*(input - options.rbfLocations(i + 1)).^2);
    end
    
elseif(strcmp(options.type,'rbf+trend'))
    %Gaussian radial basis functions + linear trend
    
    for i = 0:(options.nBasis - 3)
        f(i + 1,:) = exp(-(1/(2*options.rbfStd^2))*(input - options.rbfLocations(i + 1)).^2);
        df(i + 1,:) = -(1/options.rbfStd^2)*(input - options.rbfLocations(i + 1)).*...
            exp(-(1/(2*options.rbfStd^2))*(input - options.rbfLocations(i + 1)).^2);
    end
    %linear trend
    f(options.nBasis - 1,:) = 1;
    df(options.nBasis - 1,:) = 0;
    f(options.nBasis,:) = input;
    df(options.nBasis,:) = 1;
    
else
   error('Unknown basis function type') 
end

end

