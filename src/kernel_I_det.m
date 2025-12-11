function [f, df_dX, df_dpsi] = kernel_I_det(n, psi0, x, par)
% integral kernel and derivatives for
% deterministic contribution to valley-splitting
% INPUT: psi0 and E0 are ground-state solutions of the longitudinal confinement problem

% exponential function
  G0  = par.G0;
  arg = (n*G0(3) + 2*par.k0)*par.z;

  exp_iarg = exp(-1i*arg);

% gaussians
  arg_x   = 0.5 * n * G0(1) * par.l_x;
  arg_y   = 0.5 * n * G0(2) * par.l_y;

  arg_x_sq = arg_x*arg_x;
  arg_y_sq = arg_y*arg_y;

  gauss_x = exp( -arg_x_sq );
  gauss_y = exp( -arg_y_sq );
  gauss   = gauss_x * gauss_y;

% energy contribution from in-plane (QD) part 
  E_QD_x = 0.5 * par.const.hbar * par.omega_x * ( 0.5 - arg_x_sq);
  E_QD_y = 0.5 * par.const.hbar * par.omega_y * ( 0.5 - arg_y_sq);
  E_QD   = E_QD_x + E_QD_y;


% potential including profile variation  
  U_x = potential_modification(x, par);
  U   = par.U_F + par.U_QW + U_x;

% effective potential
  U_eff = U + E_QD;

% squared wave function (renormalize w.r.t. dz)
  psi    = psi0/sqrt(par.dz);
  psi_sq = psi .* psi;

% integral kernel
  exp_iarg_psi_sq = exp_iarg .* psi_sq;
  f   = gauss * U_eff .* exp_iarg_psi_sq;

% integral kernel derivatives  
  if nargout > 1
    df_dX   = gauss * par.dEc .* exp_iarg_psi_sq;
    df_dpsi = gauss * U_eff .* exp_iarg .* (2*psi);
  end


end

