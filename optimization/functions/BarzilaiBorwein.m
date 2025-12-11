function [x, stats, storage] = BarzilaiBorwein(x0, par, opt)
  % timing
    t_start = tic;

  % import
    print_info = opt.print_info;

  % initialize empty storage
    storage = struct;

  % counter for function evaluations
    f_eval  = 0;
    df_eval = 0;

  % print info
    if print_info >= 1
      fprintf(1,'\nStart Barzilai Borwein method ...\n')
    end
    if print_info >= 2
      fprintf(1,'iter\tf(x)\t||df(x)||\n')
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % initial step
    x  = x0;
    [f, out]  = function_f(x, par);
    [df, out] = function_df(x, out, par);
    p  = -df;

  % counter
    f_eval  = f_eval + 1;    
    df_eval = df_eval + 1;    

  % stagnation counter   
    i_stagnate = 0;  
    
  % initial step size
    alpha = opt.BB.alpha0;
%{

    abs(alpha) * scalar_product(p,p, par)
        % too large?
      if abs(alpha) * sqrt(scalar_product(p,p, par)) > par.opt.BB.dx_max
        alpha
        disp('alpha too large: reduce')
        alpha = sign(alpha) * par.opt.BB.dx_max/sqrt(scalar_product(p,p, par))

        alpha
      end
%}

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % start optimization proceduce
    iter = 0;
    continue_itertion = 1;

    while continue_itertion == 1
    % count step
      iter = iter + 1;
     
    % store old values
      x_old  = x;
      f_old  = f;
      df_old = df;

    % update
      x  = x + alpha * p;

      if sum(isnan(x))>0
        iter
        alpha
        p
        error('NaN')
      end

      if sum(isinf(x))>0
        x
        error('Inf')
      end

    % evaluate  
      [f, out]  = function_f(x, par);
      [df, out] = function_df(x, out, par);

    % counter
      f_eval  = f_eval + 1;    
      df_eval = df_eval + 1;    
  
    % search direction 
      p  = -df;

    % difference in solution and gradient 
      s = x  - x_old;
      y = df - df_old;

    % dfnorm
      dfnorm = dfnorm_function(df, par);
  
    % print info
      print_opt_progress(iter, f, dfnorm, f_old, opt)
    
    %%%%%%%%%%%%      
    % check termination conditions
      dx_update = scalar_product(s, s, par);
      df_update = scalar_product(y, y, par);
      [continue_itertion, success, i_stagnate] = check_termination_conditions(iter, f, dfnorm, dx_update, df_update, f_old, i_stagnate, opt);

    %%%%%%%%%%%%
    % user-specified inner-loop function
      [storage] = inner_loop_function(iter, alpha, x, f, df, p, out, par, opt, storage, continue_itertion);
      
    
    %%%%%%%%%%%%        
    % compute new step size
      alpha = alpha_BarzilaiBorwein(s, y, iter, opt);

    % too large?
    %{
      if abs(alpha) * sqrt(scalar_product(p,p, par)) > par.opt.BB.dx_max
        alpha
        disp('alpha too large: reduce')
        alpha = sign(alpha) * par.opt.BB.dx_max/sqrt(scalar_product(p,p, par))

        alpha
      end
    %}

  
    end
   
    % timing
      t_stop = toc(t_start);

    %%%%%%%%%%%%
    % set stats
      stats.f_eval  = f_eval;
      stats.df_eval = df_eval;
      stats.success = success;
      stats.iter    = iter;
      stats.CPU_time = t_stop;

    % print info
    if print_info > 0
      fprintf(1,'  function evaluations:   %d\n',f_eval)
      fprintf(1,'  derivative evaluations: %d\n',df_eval)
      fprintf(1,'  CPU time: %.4f sec\n',stats.CPU_time)
    end
       
   
   
      
   

end