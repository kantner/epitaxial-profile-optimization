function [S, E] = solve_schroedinger_LK(x, par)
% solve Schrödinger eigenvalue problem in longitudinal direction (growth direction)
% for prescribed epitaxial profile X = X_QW + x

  % grid
    N  = par.N;
    dz = par.dz;

  % number of eigenvalues
    neigs = par.neigs;

  % kinetic energy operator
    switch(par.discretization)
      case 1 % finite differences
        Lap = Laplacian(N);
        T   = - par.const.hbar*par.const.hbar/(2*par.mass_l) * 1/(dz^2) * Lap;
        %T_psi = @(psi)  - par.const.hbar*par.const.hbar/(2*par.mass_l) * 1/(dz^2) * ( circshift(psi,+1)-2*psi+circshift(psi,-1) );
      case 2 % spectral
        %T =   par.const.hbar*par.const.hbar/(2*par.mass_l) * par.N * ifft(ifft(diag(par.k.^2 ))')';       
        T_psi = @(psi) par.const.hbar*par.const.hbar/(2*par.mass_l) * real(ifft(par.k.^2 .* fft(psi)));
      otherwise
        error('not implemented')
    end

  % potential energy operator
    U_x = potential_modification(x, par);
    U   = par.U_QW + par.U_F + U_x;

  % find reference energy
    threshold = 0.99;    
    %idx       = find(par.QW_indicator > threshold);
    idx = (par.QW_indicator > threshold * max(par.QW_indicator));
    %if isempty(idx)
    %  U_min = 0;
    %else
      U_min = min(U(idx));
    %end
    %U_min
  %}

 % solve eigenvalue problem
    switch(par.discretization)
      case 1 % finite differences        
        H = (T + diag(sparse(U)) - U_min * speye(N))/par.energy_scale;
        
        opts.disp  = 0;
        [S, E] = eigs(H, neigs, 'smallestabs', opts);

      case 2 % spectral        
        H_psi = @(psi) (T_psi(psi) + (U - U_min).*psi)/par.energy_scale;

        opts.issym = 1; % use method for Hermitian matrix
        opts.disp  = 1;
        opts.p     = 120;
        opts.maxit = 1000;
        %opts.v0    = par.QW_indicator;    
         
        [S, E] = eigs(H_psi, N, neigs, 'smallestreal', opts);

      otherwise
        error('not implemented')
    end

  % rescale energy to SI units  
    E = diag(E) * par.energy_scale + U_min;

  % DEBUG
    DEBUG = 0;
    if DEBUG == 1
      figure(243423);clf;hold all;
      plot(par.z/par.units.nm, U/par.units.meV, 'k-')
      for i = 1 : neigs
      plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * abs(S(:,i)).^2, 'b-')
      end
      box on
      xlabel('z (nm)')
      ylabel('energy (meV)')
      title('DEBUG: Schrödinger solver')

      drawnow
      pause

    end


end

