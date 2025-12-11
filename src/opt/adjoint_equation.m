function [chi] = adjoint_equation(x, out, par)  
% solve adjoint problem using Green's function
% chi has dimension [chi] = energy^(-1) * length^(-0.5)
% INPUT
%   x ..... epitaxial profile modification
%   out ... output struct from compute_valley_splitting (which contains functional derivatives)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % RHS: sum of functional gradients w.r.t. psi
    R = zeros(par.N, 1);
    
  % contribution from J0
    if sum(abs(par.opt.a)) > 0
      %{
      DJ0_dPsi =   (out.dJ0_dM * out.dM_dnu    + out.dJ0_dV * out.dV_dnu)    * out.D_nu_dpsi ...
                 + (out.dJ0_dM * out.dM_dsigma + out.dJ0_dV * out.dV_dsigma) * out.D_sigma_dpsi;
      %}
      DJ0_dPsi =   out.dJ0_dnu    * out.D_nu_dpsi ...
                 + out.dJ0_dsigma * out.D_sigma_dpsi;
      R = R + DJ0_dPsi;
    end
  
  % contribution from J1
    if par.opt.w(1) ~= 0
      R = R + out.dJ1_dpsi;
    end

  % contribution from J2
    if par.opt.w(2) ~= 0
      R  = R + out.dJ2_dpsi;
    end

  % there is no contribution from J3  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute complete eigendecomposition of Hamiltonian
  if ismember(par.opt.adjoint_method,[1 2 3])
  % set number of eigenvalues
    par.neigs = par.N;
  
  % solve Schrödinger  
    [S, E] = solve_schroedinger(x, par);
  
  % find ground state
    psi_ref = par.QW_indicator; % set reference wave function for overlap-computation
    idx_gnd = select_ground_state(S, E, psi_ref, x, par);
  
  % check if ground state index is identical to previous computation
    if idx_gnd ~= out.idx_gnd
      error('inconsistent!')
    end
    
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate Hamiltonian
  if ismember(par.opt.adjoint_method,[4 5 6])

  % check if ground state index is identical to previous computation
    if par.discretization ~= 1
      warning('Adjoint method only implemented for finite differences. Discretization may be inconsistent')
    end
    
  % build up Hamiltonian
    % kinetic energy operator
      Lap = Laplacian(par.N);
      T   = - par.const.hbar*par.const.hbar/(2*par.mass_l) * 1/(par.dz^2) * Lap;
  
    % potential energy operator
      U_x = potential_modification(x, par);
      U   = sparse(par.U_QW + par.U_F + U_x);
  
    % Hamiltonian    
      H = T + diag(U);
    
    % LHS matrix: singular (rank-deficit) operator with wavefunction psi0 as nullspace
      A = H - out.E0 * speye(par.N);
  
  end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% built up solution using Green's function
  %fprintf(1, 'construct chi using Greens function ...')
  %solution_mode = 3; % there are three equivalent solution modes with difference performance
  % solution_mode = 3 is fast and accurate

% import
  %E       = out.E;
  %idx_gnd = out.idx_gnd;
  
  %tic
  switch(par.opt.adjoint_method)
    case 1 % iterate over eigenfunctions

      % allocate memory
        chi = zeros(par.N, 1);
    
      % import
        psi     = S/sqrt(par.dz);        
    
      % sweep over all eigenstates except for ground state
        for i = setdiff(1 : par.N, idx_gnd)    
          % compute projection of R on state
            r = par.dz * (psi(:,i)' * R);
    
          % add to chi
            chi = chi + r * psi(:,i)/(E(idx_gnd) - E(i));      
        end

    case 2 % build up Green's operator

      % index range excluding ground state
        idx_range = setdiff(1 : par.N, idx_gnd);
  
      % prepare matrices  
        S_reduced = S(:,idx_range);
        inv_diffE = diag(1./(E(idx_gnd)-E(idx_range))); % diagonal matrix


      % Green's tensor  
        G = S_reduced * inv_diffE * S_reduced';

      % matrix multiplication  
        chi = G*R;

    case 3 % evaluate Green's tensor action on RHS

      % index range excluding ground state
        idx_range = setdiff(1 : par.N, idx_gnd);
  
      % prepare matrices  
        S_reduced = S(:,idx_range);
        inv_diffE = 1./(E(idx_gnd)-E(idx_range)); % vector

      % evaluate: matrix * (vector .* (matrix * vector) )
        chi = S_reduced * ( inv_diffE .* (S_reduced'*R));
    case 4 % backslash: least-squares solution + projection
      % solve  
        chi = -A\R;

      % projection: enforce orthogonality
        chi = chi - out.psi0 * (out.psi0'*chi);

    case 5 % pseudo inverse
      % compute pseudo inverse
        pinvA = pinv(full(A));
      
      % solve  
        chi = -pinvA * R;  

    case 6 % solve augmented system
      % solve rectangular system [A; x0^T] * chi = [-R; 0] 
      % this incorporates the orthogonality condition directly
        scaling = sqrt(par.dz)^4;
        A_augmented = [A; out.psi0'*scaling];
        b_augmented = [-R; 0];

      % solve (backslash will provide a least squares solution that is
      % consistent with the Green's function solution)
        chi = A_augmented\b_augmented;
  
    otherwise
      error('not implemented')
  end
  %toc

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DEBUG = 0;
  if DEBUG == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % check if chi solves the adjoint equation
  

    switch(par.discretization)
      case 1% finite difference
      % build up Hamiltonian
      % kinetic energy operator
        Lap = Laplacian(par.N);
        T   = - par.const.hbar*par.const.hbar/(2*par.mass_l) * 1/(par.dz^2) * Lap;
    
      % potential energy operator
        U_x = potential_modification(x, par);
        U   = sparse(par.U_QW + par.U_F + U_x);
    
      % Hamiltonian    
        H = T + diag(U);
      
      % LHS matrix: singular (rank-deficit) operator with wavefunction psi0 as nullspace
        A = H - out.E0 * speye(par.N);
    
      % import
        psi     = out.S/sqrt(par.dz);  
        idx_gnd = out.idx_gnd;

      % evaluate
        LHS = A * chi;
        RHS = psi(:,idx_gnd) * (par.dz * psi(:,idx_gnd)'*R) - R;

      case 2 % spectral

      % apply Hamiltonian  
        H_chi = apply_H(chi, x, par) * par.energy_scale;

      % import
        psi     = out.S/sqrt(par.dz);        

      % evaluate  
        LHS   = H_chi - out.E0 * chi;
        RHS   = psi(:,out.idx_gnd) * (par.dz * psi(:,out.idx_gnd)'*R) - R;

      otherwise
        error('not implemented')
    end



      
  
      [LHS, RHS, LHS - RHS]
  

    % check orthogonality
      orthogonality1 = out.psi0' * chi
      orthogonality2 = out.psi0' * RHS

    %%%%%%%%%%%%  
    % plot  
      figure(123456789);clf; hold all;
        plot(par.z/par.units.nm, LHS,'ro-','DisplayName','LHS: (H-E_0)*\chi')
        plot(par.z/par.units.nm, RHS,'bo-','DisplayName','RHS: (|\psi_0><\psi_0|-I)R')
    
        legend()
        box on
        xlabel('z (nm)')
        ylabel('LHS/RHS of adjoint equation')
        drawnow
        ylabel('verify solution of adjoint equation')
        
      pause


  end

end

