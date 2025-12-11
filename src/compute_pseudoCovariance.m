function [C, dC_dX, dC_dpsi] = compute_pseudoCovariance(psi0, x, par)
% compute pseudo-covariance of complex-valued valley-coupling parameter
% and functional derivatives

    int_rand_sq    = zeros(length(par.n_range),1);
    
    if nargout > 1
      dC_dX   = zeros(par.N, 1);
      dC_dpsi = zeros(par.N, 1);
    end

  % compute for all n  
    for in = 1 : length(par.n_range)
      [f,df_dX,df_dpsi] = kernel_I_rand_sq(par.n_range(in), psi0, x, par);
      int_rand_sq(in)  = par.dz * sum(f);

      % derivatives
      if nargout > 1
        dC_dX   = dC_dX ...
                  + par.D4(in) * df_dX;
        dC_dpsi = dC_dpsi ...
                  + par.D4(in) * df_dpsi;
      end
    end     

    C  = sum(par.D4 .* int_rand_sq);


end