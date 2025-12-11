function Delta = compute_valley_splitting_det(F, x, par)
% compute deterministic components
% F is a shifted envelope belonging to the +k0 sector (column vector)


U_x   = par.dEc * x; % potential of modulation
U_tot = par.U_QW + par.U_F + U_x;



Delta = 0;
for in = 1 : length(par.n_range)

  n = par.n_range(in);

  arg_x = n * par.G0(1)*par.l_x/2;
  arg_y = n * par.G0(2)*par.l_y/2;

  arg_x_sq = arg_x*arg_x;
  arg_y_sq = arg_y*arg_y;

  gaussian_x = exp(-arg_x_sq);
  gaussian_y = exp(-arg_y_sq);
  gaussian   = gaussian_x * gaussian_y;


  U_Fm = (U_tot + 0.5*par.const.hbar*par.omega_x*(0.5-arg_x_sq) + 0.5*par.const.hbar*par.omega_y*(0.5-arg_y_sq)) .* conj(F);

  Delta = Delta ...
          + par.C2(in) * gaussian * F' * (exp(-1i*n*par.G0(3)*par.z) .* U_Fm);

end

