function [m] = compute_effective_mass_perturbation_theory(n, k, par)

% diagonalize Hamiltonian at k
  H = Hamiltonian(k, par);
  [c,E] = eig(H);
  E = diag(E) * par.energy_scale;

% allocate
  m_inv = zeros(3,3);
  
% import      
  G_list   = par.pp.G_list;
  n_states = par.pp.N_G;

% double number of states if SOI is taken into account
  if par.pp.SOI == 1
    G_list = kron(G_list, [1,1]);
    n_states = 2*n_states;
  end

% construct momentum matrix elements
%{
  p_x = zeros(n_states, n_states);
  p_y = zeros(n_states, n_states);
  p_z = zeros(n_states, n_states);

  for n = 1 : n_states
    for m = 1 : n_states 
      for l = 1 : n_states
        p_x(n,m) = c(:,m)' * par.const.hbar * G_list(1,l) * c(:,n); 
        p_y(n,m) = c(:,m)' * par.const.hbar * G_list(2,l) * c(:,n);
        p_z(n,m) = c(:,m)' * par.const.hbar * G_list(3,l) * c(:,n);
      end
    end
  end
%}

% construct momentum matrix elements p_mn for n fixed
  p = zeros(3, n_states); % array of cartesian momentum vector elements <n|p|1>, <n|p|2>, ...
  
  for m = 1 : n_states 
    for l = 1 : n_states
      p(:,m) = p(:,m) + c(l,m)' * (par.const.hbar * G_list(:,l)) * c(l,n); 
    end
  end
%}

% construct inverse effective mass tensor
  for i = 1 : 3
    for j = 1 : 3

      % diagonal term
        if i == j
          m_inv(i,j) = 1/par.const.m0;
        end

      % off-diagonal terms
        for m = 1 : n_states 
          if m ~= n % exclude n=m

            m_inv(i,j) = m_inv(i,j) ...
                         + 2/(par.const.m0^2) * (p(i,m)*p(j,m)')/(E(n)-E(m));
          end
        end
    end
  end


% compute inverse  
  m = inv(m_inv);

end