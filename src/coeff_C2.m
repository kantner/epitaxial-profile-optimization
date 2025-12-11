function C2 = coeff_C2(n_range, tol, par)
% compute bandstructure coefficients C2

  % In the bandstructure calculation, strain is taken into account without
  % explicit smallness assumption. Here it is however more convenient to
  % employ a linear approximation.

  % import
    N_G    = par.pp.N_G;
    G_list = par.pp.G_list;

  % revert to unstrained reciprocal lattice vectors via rounding
    G_list = round(G_list * par.a0/(2*pi))  * (2*pi)/par.a0; 

  %%%%%%%%%%%%%%%%%%%%%%
  % create list of G-vectors with strain to linear order
  
  % list of integer vectors (no strain)
    g_list = G_list * par.a0/(2*pi);
  
  % add linear strain
    Id = eye(3);
    for i = 1 : N_G
      g_list(:,i) = (Id - par.eps) * g_list(:,i);
    end
  
  %%%%%%%%%%%%%%%%%%%%%%
  % compute bandstructure coefficients
  
  % allocate memory
    C2 = zeros(length(n_range),1);  
  
  % scaled shift vector (times 2*pi/a0)    
    g0 = par.G0 * par.a0/(2*pi);

  % add up coefficients       
    for i1 = 1 : N_G
      g1 = g_list(:,i1);

      for i2 = 1 : N_G
        g2 = g_list(:,i2);

        % take sum (using complex conjugation of c_-(G) = c_+^*(-G)
          g = g1 + g2;
    
        % run over all integers n
          for in = 1 : length(n_range)
            n = n_range(in);
    
          % check selection rules  
            if abs(g(1) - n*g0(1) ) <= tol
              if abs(g(2) - n*g0(2) ) <= tol
                if abs(g(3) - n*g0(3) ) <= tol
                % add up coefficients  
                  C2(in) = C2(in) + par.c(i1)' * par.c(i2)';
                end
              end
            end
          end
   
      end
    end

end