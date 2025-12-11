function [M, V, dM_dnu, dM_dsigma, dV_dnu, dV_dsigma] = rician_mean_and_variance(nu, sigma)

% mean and variance
  Q  = nu/(2*sigma);
  Q2 = Q*Q;
  [f, df] = f_rician( Q2 );
  
  sqrt_pi_half = sqrt(pi/2);

  M = sqrt_pi_half * sigma * f;
  V = 2*sigma^2 + nu^2 - M^2;

% derivatives
  if nargout > 2
    dM_dnu    = sqrt_pi_half * Q * df;
    dM_dsigma = sqrt_pi_half * (f - 2 * Q2 * df);
  
    dV_dnu    = 2*(nu - M*dM_dnu);
    dV_dsigma = 2*(2*sigma - M*dM_dsigma);
  end
end