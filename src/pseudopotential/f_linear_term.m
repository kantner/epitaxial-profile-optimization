function [f] = f_linear_term(k_scaled, eps, par)
% evaluate magnitude of linear term in energy dispersion
% target is to find a k-vector where the linear term vanished to determine
% the conduction band minimum

  % generate k-vector
    k = [0 0 k_scaled]' * 2*pi/par.a0;
  
  % compute Bloch factors
    H = Hamiltonian(k, eps, par);
    [c, ~] = eig(H);
  
  % conduction band
    C = c(:,par.pp.idx_CB);
  
  % reciprocal lattice vectors
    [par] = reciprocal_lattice_vectors(eps, par);
  
    G_list = par.pp.G_list;
  
  % number of states  
    n_states = par.pp.N_G;
  
    % double everything in case of SO
    if par.pp.SOI==1
      n_states = 2*n_states;
      G_list   = kron(G_list,[1 1]);
    end
  
    
  
  % build up linear term
    scaling = 2*pi/par.a0;

    f = 0;
    for i = 1 : n_states
      f = f + real(C(i)' * (G_list(3,i) + k(3))* C(i)) * 1/scaling;
    end
    

end