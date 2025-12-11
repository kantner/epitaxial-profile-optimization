function [J0, dJ0_dnu, dJ0_dsigma] = cost_functional_J0(nu, sigma, par)
% Evaluate cost functional
%      E1 + a1*M + a2*S  +  b1*nu + b2*sigma
% J0 = -------------------------------------
%      E2 + a3*M + a4*S  +  b3*nu + b4*sigma
%
% with M = mean(E_VS)=M(nu,sigma) and
%      S = sqrt(Var(E_VS)) = sqrt(V) = S(nu,sigma)
%
% INPUT:
%   nu    ... deterministic contribution
%   sigma ... disorder-induced contribution
%

% import
  E = par.opt.E;
  a = par.opt.a;
  b = par.opt.b;

% compute mean and variance
  if nargout > 1
      [M, V, dM_dnu, dM_dsigma, dV_dnu, dV_dsigma] = rician_mean_and_variance(nu, sigma);
  else
      [M, V] = rician_mean_and_variance(nu, sigma);
  end

% cost function
  S = sqrt(V);
  denom =  E(2) + a(3)*M + a(4)*S + b(3)*nu + b(4)*sigma;
  J0    = (E(1) + a(1)*M + a(2)*S + b(1)*nu + b(2)*sigma)/denom;

% derivatives  
  if nargout > 1
    dJ0_dnu    = 1/denom * ((a(1) - J0*a(3))*dM_dnu    + (a(2) - J0*a(4))*dV_dnu/(2*S)    + (b(1) - J0*b(3)));
    dJ0_dsigma = 1/denom * ((a(1) - J0*a(3))*dM_dsigma + (a(2) - J0*a(4))*dV_dsigma/(2*S) + (b(2) - J0*b(4)));
  end
    
end
