function H = Hamiltonian(k, par)
% Generate silicon pseudopotential Hamiltonian at k assuming strain tensor eps
% following the model by Ungersboeck (2007), DOI: 10.1109/TED.2007.902880

  % import reciprocal lattice vectors
    N_G    = par.pp.N_G;
    G_list = par.pp.G_list;

  % import diamond basis vector  
    tau    = par.pp.tau;

  % precompute quantities
    kinetic_prefactor = (par.const.hbar^2)/(2*par.const.m0);

  % allocate matrix
    H = zeros(N_G, N_G);

  % allocate if SOI is active 
    if par.pp.SOI == 1
      H_SO_x = zeros(N_G, N_G);
      H_SO_y = zeros(N_G, N_G);
      H_SO_z = zeros(N_G, N_G);
    end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % assemble matrix
    for i = 1 : N_G
      
      Gi = G_list(:,i);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      % diagonal terms G_i = G_j

      %%%%%
      % kinetic energy
        H(i,i) = kinetic_prefactor * (k + Gi)' * (k + Gi) / par.energy_scale;

      %%%%%
      % local pseudopotential (diagonal)
        S      = 1; % structure factor for diagonal Gi = Gj
        H(i,i) = H(i,i) ...
                 + S * local_pseudopotential(0, par) /par.energy_scale;

      %%%%%  
      % nonlocal pseudopotential (diagonal)
        if par.pp.nonlocal == 1
          Ki      = k + Gi;
          Ki_norm = sqrt(Ki'*Ki);

          V_nloc = 4*pi/par.Omega_a * par.pp.A0 * par.pp.R0^3 * nonlocal_well_F_0_s(par.pp.R0 * Ki_norm, par.pp.R0 * Ki_norm);

          H(i,i) = H(i,i) ...
                   + S * V_nloc / par.energy_scale;
        end   

      %%%%%  
      % spin-orbit
        if par.pp.SOI == 1
        % provide B(Ki) used below  
          Bi = func_B(Ki_norm/par.pp.zeta);
        end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      % offdiagonal terms Gi =/= Gj
        for j = i+1 : N_G

          Gj     = G_list(:,j);

          Gij    = Gi - Gj;
          Gij_sq = Gij' * Gij;

        % structure factor
          S = cos(0.5 * tau' * Gij);

        %%%%%  
        % local pseudopotential  
          V_loc = local_pseudopotential(sqrt(Gij_sq), par);
          H(i,j) = S * V_loc / par.energy_scale;

        %%%%%  
        % nonlocal pseudopotential
          if par.pp.nonlocal == 1
            Kj      = k + Gj;
            Kj_norm = sqrt(Kj'*Kj);

            V_nloc = 4*pi/par.Omega_a * par.pp.A0 * par.pp.R0^3 * nonlocal_well_F_0_s(par.pp.R0 * Ki_norm, par.pp.R0 * Kj_norm);

            H(i,j) = H(i,j) ...
                     + S * V_nloc / par.energy_scale;
          end

        %%%%%%%
        % add conjugate term  (offdiagonal)
          H(j,i) = H(i,j)';

        %%%%%%% 
        % spin-orbit
          if par.pp.SOI == 1
            Bj = func_B(Kj_norm/par.pp.zeta);
    
            Ki_cross_Kj = cross(Ki, Kj);
    
            tmp = -1i*par.pp.mu * par.Omega_a^(2/3) /(pi^2) * S * Bi * Bj / par.energy_scale;
            H_SO_x(i,j) = tmp *  Ki_cross_Kj(1);
            H_SO_y(i,j) = tmp *  Ki_cross_Kj(2);
            H_SO_z(i,j) = tmp *  Ki_cross_Kj(3);
    
            H_SO_x(j,i) = H_SO_x(i,j)';
            H_SO_y(j,i) = H_SO_y(i,j)';
            H_SO_z(j,i) = H_SO_z(i,j)';
    
          end

             

        end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  % If SO-interaction is switched on, the Hamiltonian has to be tensorized
  % with a 2x2 Pauli matrix
  % Create full matrix using tensor products
    if par.pp.SOI == 1
      sigma_0 = eye(2);
      sigma_x = [0 1; 1 0];
      sigma_y = [0 -1i; 1i 0];
      sigma_z = [1 0; 0 -1];

      H = kron(H, sigma_0) ...
          + kron(H_SO_x, sigma_x) ...
          + kron(H_SO_y, sigma_y) ...
          + kron(H_SO_z, sigma_z);

    end
    



 % visualize
 %{
 tol = 1E-12;
 figure(1243);clf;hold all;
   H_sp = sparse(round(H/tol)*tol);
   figure(1244); spy(H_sp)
 %}

end