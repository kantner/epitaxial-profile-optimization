function [x, stats, storage] = gradientDescent(x0, par, opt)
  % timing
    t_start = tic;
    
  % import
    print_info = opt.print_info;

  % initialize empty storage
    storage = struct;

  % counter for function evaluations
    f_eval  = 0;
    df_eval = 0;

  % stagnation counter   
    i_stagnate = 0;
    
  % print info
    if print_info >= 1
      fprintf(1,'\nStart gradient descent method ...\n')
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % initial step  
    x  = x0;
    [f, out] = function_f(x, par);

  % counter
    f_eval = f_eval + 1;    
    
    % compute gradient    
      [df, out] = function_df(x, out, par);

    % counter
      df_eval = df_eval + 1;
     
    % set descent direction
      p = -df;
      
    % dfnorm
      dfnorm = dfnorm_function(df, par);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%
  % start optimization proceduce
    iter = 0;
    continue_itertion = 1;

    while continue_itertion == 1
    % count step
      iter = iter + 1;
  %{
    % compute gradient    
      [df, out] = function_df(x, out, par);

    % counter
      df_eval = df_eval + 1;
     
    % set descent direction
      p = -df;
      
    % dfnorm
      dfnorm = dfnorm_function(df, par);
  %}
  
    % set initial step for linesearch
      if iter == 1
      % first step  
        alpha = opt.Armijo.alpha0;
      else
      % take last accepted step size and increase by factor
        alpha = opt.Armijo.increase * alpha;
      end
    
    % store old values
      x_old  = x;
      f_old  = f;
      df_old = df;
      
    % linesearch (function evaluation at new step is done inside)
      switch(opt.gradientDescent.linesearch_mode)
        case 1
          [alpha, x, f, out, linesearch_stats] = linesearch_Armijo(alpha, x, p, f, df, par, opt);
        case 2
          [alpha, x, f, ~, out, linesearch_stats] = linesearch_Wolfe(alpha, x, p, f, df, par, opt);
        otherwise
          error('not implemented')
      end

    % counter
      f_eval  = f_eval  + linesearch_stats.f_eval;
      df_eval = df_eval + linesearch_stats.df_eval;
  
    % compute gradient    
      [df, out] = function_df(x, out, par);
      
    % counter
      df_eval = df_eval + 1;
     
    % set descent direction
      p = -df;
      
    % dfnorm
      dfnorm = dfnorm_function(df, par);
      
    % s and y
      s = x - x_old;
      y = df - df_old;

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