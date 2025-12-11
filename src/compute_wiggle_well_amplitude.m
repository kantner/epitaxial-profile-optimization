function [X_ww] = compute_wiggle_well_amplitude(X_budget, q,phi,par)
% compute amplitude of wiggle well such that Ge-budget is met

  if X_budget > 0
  % initialize
    h      = par.h_QW;
    X_ww_0 = X_budget /(0.5*( 1 + (sin(q*h-phi) +sin(phi))/(q*h) ));
  
  % set function for root finding
    mean_mode = 2;
   
    f = @(x) X_budget - mean_Ge_budget(mean_mode, x_wiggle_well(x, q, phi, par), zeros(par.N,1), par);
  
    opts = optimset('TolX',1E-36,'Display','iter');
    X_ww = fzero(f, X_ww_0,opts);
  
  else
  % set to zero directly    
    X_ww = 0;
  end

end
