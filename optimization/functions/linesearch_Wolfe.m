function [alpha, x_trial, f_trial, df_trial, out_trial, stats] = linesearch_Wolfe(alpha, x, p, f, df, par, opt)
% linesearch to satisfy (weak) Wolfe conditions
% based on Bierlaire (2015) Optimization: principles and algorithms, EPFL Press. Section 11.3
% video tutorial: https://www.youtube.com/watch?v=sXMi1D2E9QQ

% import
  maxiter     = opt.Wolfe.maxiter;
  c1          = opt.Wolfe.c1;
  c2          = opt.Wolfe.c2;
  lambda      = opt.Wolfe.lambda;
  print_info  = opt.Wolfe.print_info;

% counter
  f_eval  = 0;
  df_eval = 0;

%%%%%%%%%%%%%%%%%%%%%  
% check is input is alright
  if sum(isnan(p)) > 0
    error('p contains NaN')
  end

  if sum(isinf(p)) > 0
    error('p contains Inf')
  end

  if sum(isnan(x)) > 0
    error('x contains NaN')
  end

  if sum(isinf(x)) > 0
    error('x contains Inf')
  end

  if sum(isnan(alpha)) > 0
    error('alpha contains NaN')
  end

  if sum(isinf(alpha)) > 0
    error('alpha contains Inf')
  end  

%%%%%%%%%%%%%  
% checks 
  if 0 < c1 && c1 < c2 && c2 < 1
    % pass
  else
    % fail
      error('it must hold 0 < c1 < c2 < 1');
  end

  if lambda > 1
    % pass
  else
    % fail
      error('it must hold lambda > 1');
  end

%%%%%%%%%%%%%%%%  
% initialize counter  
  iter    = 0;

% initialize search interval
  alpha_min = 0;
  alpha_max = Inf;

% print info  
  if print_info >= 1
    fprintf(1,'\nStart Wolfe linesearch ...\n')
    fprintf(1,'iter\talpha\t\tWolfe 1\tWolfe 2\n')
  end

% evaluate scalar product for Armijo condition
  p_prod_df = scalar_product(p, df, par);  

%%%%%%%%%%%%%%%%%%%  
% linesearch
  continue_iteration = 1;
  while continue_iteration == 1
    % count iteration
      iter = iter + 1;

    % make step
      x_trial = x + alpha * p;

    % evaluate function
      [f_trial, out_trial] = function_f(x_trial, par);

    % counter  
      f_eval  = f_eval + 1;

    % print
      if print_info >= 1
        fprintf(1,'%d\t%.4e\t',iter,alpha);
      end

    % check first Wolfe condition (Armijo rule)
      if f_trial <= f + c1 * alpha * p_prod_df
        % first Wolfe condition is satisfied
        
        % print
          if print_info >= 1
            fprintf(1,'%s\t','true');
          end


        % check second Wolfe condition
          [df_trial, out_trial] = function_df(x_trial, out_trial, par);

        % counter  
          df_eval  = df_eval + 1;
  
          if -scalar_product(p, df_trial, par) <= - c2 * p_prod_df
            % second Wolfe condition is satisfied
        
            % print
              if print_info >= 1
                fprintf(1,'%s\t%s\n','true','terminate');
              end

            % terminate
              continue_iteration = 0;
              success = 1;

          else
            % second Wolfe condition is violated
            % alpha is too small
            
            % print
              if print_info >= 1
                fprintf(1,'%s\t%s\n','false','increase step');
              end


            % check for special case
              if alpha_max == Inf
                % special case where alpha_max has never been reduced to a finite value
                % increase step by factor (lambda > 1)
                  alpha = opt.Wolfe.lambda * alpha;

              else
                % regular case
                % set current alpha as new min value
                  alpha_min = alpha;

                % set new alpha in center of interval
                  alpha = 0.5*(alpha_min + alpha_max);  
              end

          end


          
        
  
      else
        % reject: first Wolfe condition is violated
        % alpha is too large

        % set current alpha as new max value
          alpha_max = alpha;

        % set new alpha in center of interval  
          alpha = 0.5*(alpha_min + alpha_max);

        % print
          if print_info >= 1
            fprintf(1,'%s\t\t%s\n','false','reduce step');
          end  

      end

    % check if maxiter is exceeded
      if iter == maxiter
        continue_iteration = 0;
        success = 0;

        if print_info >= 1
          fprintf(1,'terminate: max number of iterations reached.\n');
        end

        if print_info >= 1
        fprintf(1,'export last results (ignore that line search failed).\n');
        [f_trial, out_trial]  = function_f(x_trial, par);
        [df_trial, out_trial] = function_df(x_trial, out_trial, par);
        end
      end
   
  end

  if print_info >= 1
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