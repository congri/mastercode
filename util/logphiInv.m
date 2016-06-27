function [u] = logphiInv(X, u0)
%Inverse of logphi function,
%logphi = log(normcdf(z)).

c = [ 0.00048204; -0.00142906; 0.0013200243174; 0.0009461589032;
       -0.0045563339802; 0.00556964649138; 0.00125993961762116;
       -0.01621575378835404; 0.02629651521057465; -0.001829764677455021;
       2*(1-pi/3); (4-pi)/3; 1; 1];
   
r = [ 1.2753666447299659525; 5.019049726784267463450;
        6.1602098531096305441; 7.409740605964741794425;
        2.9788656263939928886];
    
q = [ 2.260528520767326969592;  9.3960340162350541504;
       12.048951927855129036034; 17.081440747466004316; 
        9.608965327192787870698;  3.3690752069827527677];
   
%nested function for performance
function [lp] = logphi2(z)
  z = real(z);                                 % support for real arguments only
  id1 = z.*z < 0.0492;                                 % first case: close to zero
  id2 = z < -11.3137;
  if(id1)
      lp0 = -z/sqrt(2*pi);
      f = 0; 
      for i = 1:14
         f = lp0.*(c(i) + f); 
      end
      lp = -2*f - log(2);
  elseif(id2)
      % second case: very small
      num = 0.5641895835477550741; 
      for i=1:5
          num = -z.*num/sqrt(2) + r(i);
      end
      den = 1.0;                   
      for i=1:6
          den = -z.*den/sqrt(2) + q(i);
      end
      e = num./den;
      if(norm(e) == 0)
          wanring('e == 0')
      end
      lp = log(e/2) - z.^2/2;
      
  else
      lp = log(erfc(-z/sqrt(2))/2);  % third case: rest
  end
end


u = fzero(@(u)(logphi2(u) - X), u0);


end

