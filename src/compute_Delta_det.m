function [Delta_det, D_Delta_det_dX, D_Delta_det_dpsi, Delta_det_n] = compute_Delta_det(psi0, x, par)
% compute deterministic part of complex-valued valley-splitting parameter
  
% allocate
  integral = zeros(length(par.n_range), 1);

  if nargout > 1
    D_Delta_det_dX   = zeros(par.N, 1);
    D_Delta_det_dpsi = zeros(par.N, 1);
  end
         
  for in = 1 : length(par.n_range)
    [f,df_dX,df_dpsi] = kernel_I_det(par.n_range(in), psi0, x, par);

    integral(in) = par.dz * sum(f);

    if nargout > 1
      % compute functional derivates
        D_Delta_det_dX   = D_Delta_det_dX ...
                           + par.C2(in) * df_dX;
        D_Delta_det_dpsi = D_Delta_det_dpsi ...
                           + par.C2(in) * df_dpsi;

    end
  end

% print  
%  [par.n_range', par.C2, integral/par.units.meV]
    
  Delta_det = sum(par.C2 .* integral);

% export details about contributions from resonances index by n
  if nargout >= 4    
    Delta_det_n    = par.C2 .* integral;
  end

end

