function [X_mean, dX_mean_dx, dX_mean_dpsi] = mean_Ge_budget(mode, x, psi0, par)
% compute mean Ge budget (two different approaches)

% total epitaxial profile
  par.opt.window       = par.QW_indicator;
  X = par.X_QW + par.opt.window .* x;

  switch(mode)
    case 1
      % based on wave function overlap
        X_mean = psi0' * (X .* psi0);

      % functional derivatives
        if nargout > 1
          %dX_mean_dx   = 2 * X .* psi0/sqrt(par.dz);
          %dX_mean_dpsi = (psi0.^2)/par.dz;
          dX_mean_dx   = par.opt.window .* (psi0.^2)/par.dz;
          dX_mean_dpsi = 2 * X .* psi0/sqrt(par.dz);
        end

    case 2
      % based on QW indicator
        X_mean = sum(par.QW_indicator .* X)/sum(par.QW_indicator);

      % functional derivatives
        if nargout > 1
          dX_mean_dx   = par.QW_indicator/(par.dz * sum(par.QW_indicator) );
          dX_mean_dpsi = zeros(par.N, 1);
        end

    otherwise
      error('not implemented')
  end

end