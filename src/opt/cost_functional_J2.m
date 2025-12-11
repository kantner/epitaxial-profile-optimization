function [J2, dJ2_dx, dJ2_dpsi] = cost_functional_J2(x, psi0, par)
% cost functional enforcing target Ge concentration

  if par.opt.w(2) ~= 0
    
    % compute mean in mode 2
      mode = 2;
      if nargout == 1
        [X_mean] = mean_Ge_budget(mode, x, psi0, par);
      else
        [X_mean, dX_mean_dx, dX_mean_dpsi] = mean_Ge_budget(mode, x, psi0, par);
      end
      
    % cost functional
      J2 = 0.5 * par.opt.w(2) * (X_mean - par.opt.X_budget)^2;
    
    % gradients
      if nargout > 1
        dJ2_dx   = par.opt.w(2) * (X_mean - par.opt.X_budget) .* dX_mean_dx;
        dJ2_dpsi = par.opt.w(2) * (X_mean - par.opt.X_budget) .* dX_mean_dpsi;
      end
     
  else 
    % if w(2) = 0 everything is zero
    J2 = 0;
    if nargout > 1
      dJ2_dx   = zeros(par.N, 1);
      dJ2_dpsi = zeros(par.N, 1);
    end
  end

end