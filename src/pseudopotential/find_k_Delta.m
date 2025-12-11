function [k_Delta] = find_k_Delta(par)
% two available modes
%   mode 1: find root of linear term in energy dispersion
%   mode 2: find minimum of energy band

% select mode
  mode = 2;

  switch(mode)
    case 1
    % find root of linear term

    % initial guess
      k_scaled_init = 0.85;
      
      opts = optimset('TolX',1E-16,'Display','iter');
      k_scaled_min = fzero(@(k_scaled) f_linear_term(k_scaled, par), k_scaled_init, opts);
    
      k_Delta = [0 0 k_scaled_min]' * 2*pi/par.a0;
    
      %k0 = k_Delta(3)*par.a0/(2*pi);

    case 2
    % scalar minimization of energy band using interpolation extrapolation

      plot_progress = 0; % 0 = off | 1 = on

      tol  = 1E-7;
      nmax = 1000; % max iterations

      fprintf('\n')
      fprintf('compute k_Delta ...\n')
      fprintf('iter\tk_0 (2*pi/a0)\tE/E_scale\tres\n')

      k_scaled = zeros(3, 1);
      E_scaled = zeros(3, 1);

    % select initial values bracketing the minimum
      k_scaled(1) = 0.1;
      k_scaled(2) = 0.9;

    % choose midpoint as third value
      k_scaled(3) = 0.5*(k_scaled(1) + k_scaled(2));

    % evaluate energies
      for i = 1 : 3
      k = [0 0 k_scaled(i)]' * 2*pi/par.a0;
      H = Hamiltonian(k, par);
      E = eig(H);
      E_scaled(i) = E(par.pp.idx_CB);
      end
      

      if plot_progress == 1
        figure(200); clf; hold all;
          plot(k_scaled,E_scaled,'ko')
          box on
          xlabel('k')
          ylabel('E')
          drawnow
      end


      continue_iteration = 1;
      iter = 0;

      while continue_iteration == 1
      % count iteration  
        iter = iter + 1;

      % quadratic fit E(k) = c1 k^2 + c2 k + c3
        A = zeros(3,3);
        for i = 1 : 3
        A(i,:) = [k_scaled(i)^2 k_scaled(i) 1];
        end
        c = A\E_scaled;

        if plot_progress == 1
          k_range = linspace(min(k_scaled), max(k_scaled), 51);
        
          figure(200); hold all;
            plot(k_range,c(1) * k_range.^2 + c(2) * k_range + c(3),'r-')
            drawnow
        end

      % store previous update to check convergence
        if iter > 1
          k_scaled_new_prev = k_scaled_new;
        end
        
      % interpolate: new minimum candidate      
        k_scaled_new = - c(2)/(2*c(1));

      % energy at new minimum candidate
        k = [0 0 k_scaled_new]' * 2*pi/par.a0;
        H = Hamiltonian(k, par);
        E = eig(H);
        E_scaled_new = E(par.pp.idx_CB);
  
        if plot_progress == 1        
          figure(200); hold all;
            plot(k_scaled_new,E_scaled_new,'gx')
            drawnow
        end  


      % select two lowest energy points from previous run and new point
        [~,ord] = sort(E_scaled,'ascend');
  
        k_scaled = [k_scaled(ord(1:2)); k_scaled_new];
        E_scaled = [E_scaled(ord(1:2)); E_scaled_new];
  

      % residuum and print
        if iter > 1
        % choose residuum
          res_mode = 2;
          switch(res_mode)
            case 1 % difference of k-vectors
              res = abs(k_scaled_new - k_scaled_new_prev);
            case 2 % width of k window
              res = max(k_scaled) - min(k_scaled);
          end

          fprintf('%d\t%.6f\t%.6f\t%.4e\n', iter, k_scaled_new, E_scaled_new, res);
        else
          fprintf('%d\t%.6f\t%.6f\n', iter, k_scaled_new, E_scaled_new);
        end


      % check termination conditions  
        if iter == nmax
          message = 'stop: maximum number of iterations reached.';
          continue_iteration = 0;
          k0_scaled = k_scaled_new;
        end

        %{
          if max(E_scaled) - min(E_scaled) < tol
          message= 'stop: energy differences below tol.';
          continue_iteration = 0;
          k0_scaled = k_scaled_new;
        end
        %}
        
        if iter > 1
        if abs(k_scaled_new_prev - k_scaled_new) < tol
          message = 'stop: k vector update below tol.';
          continue_iteration = 0;
          k0_scaled = k_scaled_new;
        end
        end

        if continue_iteration == 0
          disp(message)
        end

      end

    % export
      k_Delta = [0 0 k0_scaled]' * 2*pi/par.a0;



  end

end