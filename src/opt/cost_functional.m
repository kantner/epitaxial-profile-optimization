function [out] = cost_functional(x, compute_derivatives, par)

  %%%%%%%%%%%%%%%%%%%%%
  % apply filter and window
  % Replace x --> window * conv(filter, x) at the beginning of the
  % functional computation. This modification is passed to all subroutines.


  %%%%%%%%%%%%%%%%%%%%%
  % check if x is alright
    if sum(isnan(x)) > 0
      error('x contains NaN (before filtering)')
    end
    if sum(isinf(x)) > 0
      error('x contains Inf (before filtering)')
    end

  %%%%%%%%%%%%%%%%%%%%%
  % apply filter and window
    if isfield(par,'opt')
    % apply filter
      if par.opt.apply_filter == 1
        [x] = apply_filter(x, par.opt.filter);
      end
  
    % apply window
      if par.opt.apply_window == 1
        [x] = apply_window(x, par.opt.window);
      end
    end

  %%%%%%%%%%%%%%%%%%%%%
  % check if x is alright
    if sum(isnan(x)) > 0
      error('x contains NaN (after filtering)')
    end
    if sum(isinf(x)) > 0
      error('x contains Inf (after filtering)')
    end

  %%%%%%%%%%%%%%%%%%%%%
  % compute valley splitting
    [out] = compute_valley_splitting(x, compute_derivatives, par);

  % import
    %M     = out.M; % mean(E_VS)
    %V     = out.V; % var(E_VS)
    nu    = out.nu;
    sigma = out.sigma;
    psi0  = out.psi0;

  %%%%%%%%%%%%%%%%%%%%%  
  % J0: figure of merit    
    if compute_derivatives == 1
      [J0, dJ0_dnu, dJ0_dsigma] = cost_functional_J0(nu, sigma, par);
    else
      [J0] = cost_functional_J0(nu, sigma, par);
    end

  %%%%%%%%%%%%%%%%%%%%%
  % J1: Ge budget (1)
    if compute_derivatives == 1
      [J1, dJ1_dx, dJ1_dpsi] = cost_functional_J1(x, psi0, par);
    else
      [J1] = cost_functional_J1(x, psi0, par);
    end

  %%%%%%%%%%%%%%%%%%%%%
  % J2: Ge budget (2)
    if compute_derivatives == 1
      [J2, dJ2_dx, dJ2_dpsi] = cost_functional_J2(x, psi0, par);
    else
      [J2] = cost_functional_J2(x, psi0, par);
    end

  %%%%%%%%%%%%%%%%%%%%%
  % J3: admissible range
    if compute_derivatives == 1
      [J3, dJ3_dx] = cost_functional_J3(x, par);
    else
      [J3] = cost_functional_J3(x, par);
    end
    
  %%%%%%%%%%%%%%%%%%%%%
  % total cost
    J_tot = J0 + J1 + J2 + J3;
  
  %%%%%%%%%%%%%%%%%%%%%  
  % compile output structure (add to out-struct from valley splitting computation above)
    out.J_tot = J_tot;
    out.J0    = J0;
    out.J1    = J1;
    out.J2    = J2;
    out.J3    = J3;

    
  % store derivatives for gradient computation in out-struct  
    if compute_derivatives == 1
      %out.dJ0_dM = dJ0_dM;
      %out.dJ0_dV = dJ0_dV;
      out.dJ0_dnu    = dJ0_dnu;
      out.dJ0_dsigma = dJ0_dsigma;

      out.dJ1_dx = dJ1_dx;
      out.dJ2_dx = dJ2_dx;
      out.dJ3_dx = dJ3_dx;

      out.dJ1_dpsi = dJ1_dpsi;
      out.dJ2_dpsi = dJ2_dpsi; % this is identical zero ...
    end

    % the following items are already contained in out
    %out.nu    = nu;
    %out.sigma = sigma;
    %out.psi0  = psi0;
    %out.E0    = E0;
    %out.S       = S;
    %out.E       = E;
    %out.idx_gnd = idx_gnd;


end