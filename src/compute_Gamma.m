function [Gamma, dGamma_dX, dGamma_dpsi] = compute_Gamma(psi0, x, par)
% compute covariance of complex-valued valley-splitting parameter
% and functional derivatives

    int_rand_abs    = zeros(length(par.n_range),1);
    
    if nargout > 1
      dGamma_dX   = zeros(par.N, 1);
      dGamma_dpsi = zeros(par.N, 1);
    end

  % compute for all n  
    for in = 1 : length(par.n_range)
      [f,df_dX,df_dpsi] = kernel_I_rand(par.n_range(in), psi0, x, par);
      int_rand_abs(in)  = par.dz * sum(f);

      % derivatives
      if nargout > 1
        dGamma_dX   = dGamma_dX ...
                      + par.C4(in) * df_dX;
        dGamma_dpsi = dGamma_dpsi ...
                      + par.C4(in) * df_dpsi;
      end
    end     

    Gamma  = sum(par.C4 .* int_rand_abs);

    %{
     [par.n_range, par.C4m, int_rand_abs/par.units.meV^2 ]
     abs_Delta_rand_sq/par.units.meV^2
      pause
    %}

    if Gamma < 0
      warning(['Gamma is negative: ',num2str(Gamma,'%.4e')])
      Gamma = par.Gamma_failure;
      warning(['----> replace by dummy value: ',num2str(Gamma,'%.4e')])
    end


  %%%%%%%%%%%  
  % everthing should be real (filter out residual imaginary parts from inexact arithmetics)
    Gamma = real(Gamma);
    
    if nargout > 1
      dGamma_dX   = real(dGamma_dX);
      dGamma_dpsi = real(dGamma_dpsi);
    end




end