function [x, stats, storage] = LBFGS(x0, par, opt)
  % timing
    t_start = tic;

  % import
    print_info = opt.print_info;
    m          = opt.LBFGS.memory; % L-BFGS memory depth

  % initialize empty storage
    storage = struct;

  % counter for function evaluations
    f_eval  = 0;
    df_eval = 0;

  % print info
    if print_info >= 1
      fprintf(1,'\nStart L-BFGS method (m=%d) with Wolfe linesearch ...\n',m)
    end

  % allocate memory for history storage (this has nothing to do with the storage struct)
    s = zeros(par.N, m); % store solution updates s = x(k) - x(k-1)
    y = zeros(par.N, m); % store gradient updates y = df(x(k)) - df(x(k-1))

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % initial step  
    x  = x0;
    [f, out]  = function_f(x, par);
    %f  = function_f(x, par);

  % counter 
    f_eval = f_eval + 1;

  % stagnation counter   
    i_stagnate = 0;
    
  % hist size counter
    hist_size = 0;

  % skip counter
    skip_counter = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % start optimization proceduce
    iter = 0;
    continue_itertion = 1;

    while continue_itertion == 1
    % count step
      iter = iter + 1;

    % compute search direction p
      if iter == 1  
      % The first step is different, as there is no history of gradients yet.
      % Therefore, we choose as simple gradient descend step initially.

      % compute gradient    
        %df = function_df(x, par);
        [df, out] = function_df(x, out, par);
  
      % counter   
        df_eval = df_eval + 1;
  
      % dfnorm
        dfnorm = dfnorm_function(df, par);  
  
      % simple gradient descend in first step
        p = - df;
  
      else % iter > 1
      % evaluate p = - B*df using history of s and y and two-loop-recursion

        % set range of history-arrays
        %{
          if iter <= m
          % arrays are still filled up  
            range = 1 : iter-1;
          elseif iter > m
          % arrays are full  
            range = 1 : m;
          end
        %}
        range = 1 : hist_size;
        
        % evaluate p  
          p = - two_loop_recursion(df, s(:,range), y(:,range));
      end
          
    % set initial step for linesearch
      if iter == 1
        alpha = opt.Wolfe.alpha_init;
      else
        alpha = opt.Wolfe.alpha0;
      end

      if iter > 1
        if alpha ~= 1.0
          warning('alpha should be set to 1 for L-BFGS')
          fprintf(1,'set alpha to 1');
          alpha = 1.0;
        end
      end
      
    % store old values
      x_old  = x;
      f_old  = f;
      df_old = df;

    % linesearch (function evaluation at new step is done inside)
      [alpha, x, f, df, out, linesearch_stats] = linesearch_Wolfe(alpha, x, p, f, df, par, opt);

    % counter
      f_eval  = f_eval  + linesearch_stats.f_eval;
      df_eval = df_eval + linesearch_stats.df_eval;

    % dfnorm
      dfnorm = dfnorm_function(df, par);  
  
    %%%%%%%%%%%%%%%%%  
    % history storage
      s_new = x - x_old;
      y_new = df - df_old;

      if s_new' * y_new >= par.opt.LBFGS.min_curv
      % reset skip counter
        skip_counter = 0;

      % add to history  
        if hist_size < m
          hist_size = hist_size + 1;
          % fill up arrays
          % add new columns
          s(:,hist_size) = s_new;
          y(:,hist_size) = y_new;
        else

          % swap memory
          %{
          % this seems to be slow
            s(:,1:m-1) = s(:,2:m);
            y(:,1:m-1) = y(:,2:m);
          %}
          % this seems to be more efficient
          %
          s = circshift(s,[0,-1]);
          y = circshift(y,[0,-1]);
          %}

          % add new entries as last column
          s(:,m) = s_new;
          y(:,m) = y_new;
        end

      else  
          fprintf('skip L-BFGS update because curvature condition if degenerate sy=0\n')
          skip_counter = skip_counter + 1;
      end


      

    % print info
      print_opt_progress(iter, f, dfnorm, f_old, opt)
  
    %%%%%%%%%%%%      
    % check termination conditions
      dx_update = scalar_product(s_new, s_new, par);
      df_update = scalar_product(y_new, y_new, par);

      [continue_itertion, success, i_stagnate] = check_termination_conditions(iter, f, dfnorm, dx_update, df_update, f_old, i_stagnate, opt);

      if continue_itertion == 1
        if skip_counter >= par.opt.LBFGS.max_skip
        % maximum skips of history update reache
          continue_itertion = 0;
          termination = 1;
          success     = 0;
          message     = 'Reached maximum number of history update skips in a row.';
  
          fprintf(1,'Optimization has stopped after %d iterations: ', iter)
          fprintf(1,'%s\n',message);
          pause(1)
        end
      end



    %%%%%%%%%%%%
    % user-specified inner-loop function
      [storage] = inner_loop_function(iter, alpha, x, f, df, p, out, par, opt, storage, continue_itertion);

    end

    % timing
      t_stop = toc(t_start);

    %%%%%%%%%%%%
    % set stats
      stats.f_eval   = f_eval;
      stats.df_eval  = df_eval;
      stats.success  = success;
      stats.iter     = iter;
      stats.CPU_time = t_stop;

    % print info
    if print_info > 0
      fprintf(1,'  function evaluations:   %d\n',stats.f_eval)
      fprintf(1,'  derivative evaluations: %d\n',stats.df_eval)
      fprintf(1,'  CPU time: %.4f sec\n',stats.CPU_time)
    end
   
   
   

end