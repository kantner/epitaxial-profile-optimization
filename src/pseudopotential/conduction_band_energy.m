function Ec = conduction_band_energy(k_scaled, par)

  k = [0 0 k_scaled]' * 2*pi/par.a0;

  H = Hamiltonian(k, par);
  [~,E] = eig(H);
  E = diag(E);
  Ec = E(par.pp.idx_CB);

end