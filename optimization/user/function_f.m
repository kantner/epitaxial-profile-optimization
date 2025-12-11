function [f, out] = function_f(x, par)
  %model = par.model;
%{
  switch(model)
    case 1
      f = f_rosenbrock(x, par);
    case 2
      f = f_styblinksi_tang(x, par);      
    case 3
%}
      compute_derivatives = 1;  
      [out] = cost_functional(x, compute_derivatives, par);  
      f = out.J_tot;
    
%{
otherwise
      error('not implemented')
  end
%}

  if isnan(f)
    error('f evaluated to NaN')
  end

  if isinf(f)
    error('f evaluated to Inf')
  end

end