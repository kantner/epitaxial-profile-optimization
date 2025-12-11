function [tau] = strained_tau(par)

  Id = eye(3);

  tau = 0.25 * par.a0 * [1 1 1]';

%%%%%%%%%%%%%%%%
% macroscopic strain: barycenter
  tau_0 = (Id + par.eps) * tau;

%%%%%%%%%%%%%%%%  
% internal strain: circumcenter (equal distance to neighboring atoms of other fcc-sublattice)

% import strained primitive lattice vectors
  for i = 1 : 3
    a{i} = par.pp.a{i};
  end

  A = [a{1}'; a{2}'; a{3}'];
  b = zeros(3,1);
  for i = 1 : 3
    b(i) = 0.5 * a{i}'*a{i};
  end

  tau_1 = A\b;

%%%%%%%%%%%%%%%%
% convex combination
  zeta = par.pp.internal_strain;
  tau  = zeta * tau_1 + (1-zeta) * tau_0;

end