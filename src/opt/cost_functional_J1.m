function [J1, dJ1_dx, dJ1_dpsi] = cost_functional_J1(x, psi0, par)
% cost functional enforcing target Ge concentration

  if par.opt.w(1) ~= 0
    
    % compute mean in mode 1
      mode = 1;
      if nargout == 1
        [X_mean] = mean_Ge_budget(mode, x, psi0, par);
      else
        [X_mean, dX_mean_dx, dX_mean_dpsi] = mean_Ge_budget(mode, x, psi0, par);
      end
      
    % cost functional
      J1 = 0.5 * par.opt.w(1) * (X_mean - par.opt.X_budget)^2;
    
    % gradients
      if nargout > 1
        dJ1_dx   = par.opt.w(1) * (X_mean - par.opt.X_budget) .* dX_mean_dx;
        dJ1_dpsi = par.opt.w(1) * (X_mean - par.opt.X_budget) .* dX_mean_dpsi;
      end
     
  else 
    % if w(1) = 0 everything is zero
    J1 = 0;
    if nargout > 1
      dJ1_dx   = zeros(par.N, 1);
      dJ1_dpsi = zeros(par.N, 1);
    end
  end

end