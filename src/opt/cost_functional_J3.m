function [J3, dJ3_dx] = cost_functional_J3(x, par)
% cost functional enforcing Ge concentration in admissible range x_min <= x <= x_max

  if par.opt.w(3) ~= 0
    
    % total profile
      X = par.X_QW + x;
    
    % evaluate penalty function
      if nargout == 1
        [f] = f_penalty(X, par);
      else
        [f, df] = f_penalty(X, par);
      end

    % cost functional
      J3 = par.opt.w(3) * sum(f) * par.dz/par.h_QW ;
    
    % gradients
      if nargout > 1
        dJ3_dx   = par.opt.w(3) * df/par.h_QW;
      end
     
  else 
    % if w(3) = 0 everything is zero
    J3 = 0;
    if nargout > 1
      dJ3_dx   = zeros(par.N, 1);
    end
  end

end