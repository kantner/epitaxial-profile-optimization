function [alpha, x_trial, f_trial, out_trial, stats] = linesearch_Armijo(alpha, x, p, f, df, par, opt)
% linesearch based on Armijo rule
% the function

% import
  maxiter     = opt.Armijo.maxiter;
  c1          = opt.Armijo.c1;
  shrink      = opt.Armijo.shrink;
  print_info  = opt.Armijo.print_info;

% initialize counter  
  iter    = 0;

% counter
  f_eval  = 0;
  df_eval = 0;

% print info  
  if print_info == 1
    fprintf(1,'\nStart Armijo linesearch ...\n')
    fprintf(1,'iter\talpha\n')
  end

% evaluate scalar product for Armijo condition
  p_prod_df = scalar_product(p, df, par);
  

% plot linesearch  
  if opt.Armijo.plot_linesearch  == 1
    linesearch_fig = 999;  
    alpha_range = par.opt.Armijo.plot_alpha_range;

    phi_range = zeros(size(alpha_range));
    for i = 1 : length(alpha_range)
      phi_range(i) = function_f(x + alpha_range(i) * p, par);
    end
    
    figure(linesearch_fig); clf; hold all;
      xlabel('\alpha')
      ylabel('\phi(\alpha)')
      title('Armijo linesearch')
      box on
      legend()
      ylim([-1 1]*100)

      yline(f, 'g--', 'DisplayName', '\phi(0)','LineWidth',2)
      plot(alpha_range, phi_range, 'ko-', 'DisplayName', '\phi(\alpha)','LineWidth',2)
      plot(alpha_range, f + alpha_range * p_prod_df, 'r.-', 'DisplayName', '\phi(0) + \alpha \phi\prime(0)','LineWidth',2)

      set(gca,'XScale',par.opt.Armijo.plot_alpha_scale)

      drawnow


  end


%%%%%%%%%%%%%%%%%%%  
% linesearch
  continue_iteration = 1;
  while continue_iteration == 1
    % count iteration
      iter = iter + 1;

    % make step
      x_trial = x + alpha * p;

    % evaluate  
      [f_trial, out_trial] = function_f(x_trial, par);  

    % counter
      f_eval  = f_eval + 1;  

    % plot  
      if opt.Armijo.plot_linesearch  == 1
        figure(linesearch_fig); hold all;
          plot(alpha, f_trial, 'rs', 'DisplayName', 'new point')
      end

    % print
      if print_info == 1
        fprintf(1,'%d\t%.4g\t',iter,alpha);
      end

    % check Armijo rule  
      if f_trial <= f + c1 * alpha * p_prod_df
        % accept
        
        % terminate
          continue_iteration = 0;
          success = 1;

        % print
          if print_info == 1
            fprintf(1,' ... %s\n','success');
          end
  
      else
        % reject: reduce step size
          alpha = shrink * alpha;

        % print
          if print_info == 1
            fprintf(1,' ... %s\n','failed: reduce step');
          end  
      end

    % check if maxiter is exceeded
      if iter == maxiter
        continue_iteration = 0;
        success = 0;
      end
   
  end



% plot tangent of final point  
  if opt.Armijo.plot_linesearch  == 1

    % evaluate derivative (only for plotting)
      [df_trial, ~] = function_df(x_trial, out_trial, par);  
    
    figure(linesearch_fig); hold all;
    
      phi_prime = par.dz*(p'*df_trial);

      plot(alpha_range, f_trial + phi_prime * (alpha_range-alpha), 'm-', 'DisplayName', 'tangent','LineWidth',2)
      drawnow
      pause


  end



  if print_info > 0
    fprintf(1,'function evaluations:   %d \n',f_eval);
    fprintf(1,'derivative evaluations: %d \n',df_eval);
  end

%%%%%%%%%%%%%%%%%%%  
% compile stats  
  stats.iter    = iter;
  stats.success = success;
  stats.f_eval  = f_eval;
  stats.df_eval = df_eval;  

end