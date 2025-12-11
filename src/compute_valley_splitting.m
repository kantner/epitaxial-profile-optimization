function [out] = compute_valley_splitting(x, compute_derivatives, par)
% compute valley splitting for given epitaxial profile
% INPUT
%   x ... epitaxial profile modification
%   compute_derivatives = 0 (false) or 1 (true)

  %%%%%%%%%%%%%%%%%%%%
  % solve Schrödinger eigenvalue problem
    [S, E] = solve_schroedinger(x, par);

  %%%%%%%%%%%%%%%%%%%%
  % select ground state
    psi_ref = par.QW_indicator; % set reference wave function for overlap-computation
    idx_gnd = select_ground_state(S, E, psi_ref, x, par);
    
  % extract ground state wave function and  energy
    psi0 = S(:,idx_gnd);
    E0   = E(idx_gnd);

  %%%%%%%%%%%%%%%%%%%%%%%  
  % plot / debug
    plot_solution = 0;
    if plot_solution == 1
      figure(1344); clf; hold all;
        U_tot = par.U_QW + par.U_F + potential_modification(x, par);
        plot(par.z/par.units.nm, U_tot/par.units.meV, 'k-')
        for i = 1 : par.neigs   
          %if E(i) <= max_U
          plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * abs(S(:,i)).^2, '-','Color',[1 1 1]*0.5,'LineWidth',1)
          %end
        end
        i = idx_gnd;
        plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * abs(S(:,i)).^2, '-','Color',[1 0 0],'LineWidth',2)
    
        drawnow
        pause

    end

    
  %%%%%%%%%%%%%%%%%%%%%%%%  
  % compute valley splitting
    
  % determinsitic part
    if compute_derivatives == 1
      [Delta_det, D_Delta_det_dX, D_Delta_det_dpsi, Delta_det_n] = compute_Delta_det(psi0, x, par);
    else
      [Delta_det, ~, ~, Delta_det_n] = compute_Delta_det(psi0, x, par);
    end

  % random part: covariance
    if compute_derivatives == 1
      [Gamma,  dGamma_dX,  dGamma_dpsi]  = compute_Gamma(psi0, x, par);
    else
      [Gamma]  = compute_Gamma(psi0, x, par);
    end

  % random part: pseudo-covariance
    if compute_derivatives == 1
      [C,  dC_dX,  dC_dpsi]  = compute_pseudoCovariance(psi0, x, par);
    else
      [C]  = compute_pseudoCovariance(psi0, x, par);
    end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot
  %{
    pdf_mode=2; % 0 = none | 1 = circular | 2 = complex normal
    plot_complex_Delta(Delta_det, Gamma, C, pdf_mode, par)
  %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  %%%%%%%%%%%%%%%%%%%%%%%%  
  % Rice distribution parameters
    nu    = 2 * abs(Delta_det);
    sigma = sqrt(2 * Gamma);
  
  % consistency check
    if or(sigma < 0, imag(sigma)>0)
      error('sigma must be real and non-negative.')
    end

  % provide functional derivatives
    if compute_derivatives == 1
      D_nu_dpsi = 2 * real( Delta_det'/abs(Delta_det) * D_Delta_det_dpsi );
      D_nu_dX   = 2 * real( Delta_det'/abs(Delta_det) * D_Delta_det_dX    );

      D_sigma_dpsi = 1/sigma * dGamma_dpsi;
      D_sigma_dX   = 1/sigma * dGamma_dX;
    end


  %%%%%%%%%%%%%%%%%%%%%%%%
  % compute mean and variance of rice distribution  
    if compute_derivatives == 1
      [M, V, dM_dnu, dM_dsigma, dV_dnu, dV_dsigma] = rician_mean_and_variance(nu, sigma);
    else
      [M, V] = rician_mean_and_variance(nu, sigma);
    end
    % here M = mean(E_VS) and V = var(E_VS)
        
  %%%%%%%%%%%%%%%%%%%%%%%%
  % compile output
    out.M     = M;
    out.V     = V;
    
    out.Delta_det = Delta_det;
    out.Gamma = Gamma;
    out.C     = C;
    
    out.nu    = nu;
    out.sigma = sigma;
    

    out.Delta_det_n = Delta_det_n;
    
    out.psi0 = psi0;
    out.E0   = E0;

    out.S       = S;
    out.E       = E;
    out.idx_gnd = idx_gnd;

    %out.real_Delta_rand_sq = real_Delta_rand_sq; % optional

    if compute_derivatives == 1
      out.dM_dnu    = dM_dnu;
      out.dM_dsigma = dM_dsigma;
      out.dV_dnu    = dV_dnu;
      out.dV_dsigma = dV_dsigma;
  
      out.D_nu_dpsi    = D_nu_dpsi;
      out.D_nu_dX      = D_nu_dX;

      out.D_sigma_dpsi = D_sigma_dpsi;
      out.D_sigma_dX   = D_sigma_dX;

      %out.D_Delta_det_dX          = D_Delta_det_dX;
      %out.D_Delta_det_dpsi        = D_Delta_det_dpsi;
  
      %out.dabs_Delta_rand_sq_dX   = dabs_Delta_rand_sq_dX;
      %out.dabs_Delta_rand_sq_dpsi = dabs_Delta_rand_sq_dpsi;
    end

end    