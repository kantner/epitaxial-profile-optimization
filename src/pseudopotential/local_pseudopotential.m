function [V] = local_pseudopotential(q, par)
% check is non-negative
  if q < 0
    error('q must be non-negative')
  end

% import
  kF = par.pp.kF;

% rescale
  q_scale = 2*pi/par.a0;
  q       = q/q_scale;  

  kF      = kF/q_scale;  

% import weights
  w = par.pp.w;

  if 0 <= q && q < sqrt(3)
    V = ((w(1)*q + w(2))*q + w(3))*q + w(4);
  elseif sqrt(3) <= q && q < sqrt(8)
    V = ((w(5)*q + w(6))*q + w(7))*q + w(8);
  elseif sqrt(8) <= q && q < sqrt(11)
    V = ((w(9)*q + w(10))*q + w(11))*q + w(12);
  elseif sqrt(11) <= q && q < 3*kF
    V = ((w(13)*q + w(14))*q + w(15))*q + w(16);
  elseif q >= 3*kF
    V = 0;
  end

% rescale
  V = V * par.units.Ry;

end