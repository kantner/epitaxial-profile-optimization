function [effective_mass] = compute_effective_mass_finite_difference(n, k, dk, par)
% compute hessian from band dispersion E_n(k) using finite difference approximation
% n  ... band index
% k  ... k-vector
% dk ... finite difference

  
 %% energy at k
    H = Hamiltonian(k, par);

    [~,E] = eig(H);
    E = diag(E)*par.energy_scale;

    En = E(n);

 %% gradient at k: central finite difference
  %{
    gradient = zeros(3,1);

    for i = 1 : 3

    % unit vector  
      ei = zeros(3,1);
      ei(i) = 1; 

    % pertubation +dk/2 * ei
      H = Hamiltonian(k + 0.5*dk*ei, eps, par);
      [~,E] = eig(H);
      E = diag(E)*par.energy_scale;
      Ep = E(n);

    % pertubation -dk/2 * ei
      H = Hamiltonian(k - 0.5*dk*ei, eps, par);
      [~,E] = eig(H);
      E = diag(E)*par.energy_scale;
      Em = E(n);

    % gradient
      gradient(i) = (Ep-Em)/dk;

    end

    gradient

    k' * gradient 
  %}

   
 %% hessian at k
    hessian  = zeros(3,3);

    for i = 1 : 3
      % unit vector  
        ei = zeros(3,1);
        ei(i) = 1; 

      for j = 1 : 3

        % unit vector  
          ej = zeros(3,1);
          ej(j) = 1; 
    
        % pertubation + dk/2 * ei + dk/2 * ej
          H = Hamiltonian(k + 0.5*dk*ei + 0.5*dk*ej, par);
          [~,E] = eig(H);
          E = diag(E)*par.energy_scale;
          E1 = E(n);

        % pertubation + dk/2 * ei - dk/2 * ej
          H = Hamiltonian(k + 0.5*dk*ei - 0.5*dk*ej, par);
          [~,E] = eig(H);
          E = diag(E)*par.energy_scale;
          E2 = E(n);

        % pertubation - dk/2 * ei + dk/2 * ej
          H = Hamiltonian(k - 0.5*dk*ei + 0.5*dk*ej, par);
          [~,E] = eig(H);
          E = diag(E)*par.energy_scale;
          E3 = E(n);

        % pertubation - dk/2 * ei - dk/2 * ej
          H = Hamiltonian(k - 0.5*dk*ei - 0.5*dk*ej, par);
          [~,E] = eig(H);
          E = diag(E)*par.energy_scale;
          E4 = E(n);          

        % hessian
          hessian(i,j) = (E1-E2-E3+E4)/(dk*dk);

      end
          
      
    end

 %% effective mass tensor
    inv_hessian = eye(3)/hessian;
    effective_mass = par.const.hbar^2 * inv_hessian;


end

