function [S, E] = solve_schroedinger(x, par)
% solve Schrödinger eigenvalue problem in longitudinal direction (growth direction)
% for prescribed epitaxial profile X = X_QW + x

  switch(par.envelope_model)
    case 1 % Luttinger-Kohn type
      [S, E] = solve_schroedinger_LK(x, par);
    case 2 % Burt-Foreman type
      [S, E] = solve_schroedinger_BF(x, par);
    otherwise
      error('invalid option')
  end

end

