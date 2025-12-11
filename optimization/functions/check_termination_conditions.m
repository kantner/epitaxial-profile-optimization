function [continue_itertion, success, i_stagnate] = check_termination_conditions(iter, f, dfnorm, dx_update, df_update, f_old, i_stagnate, opt)

  % initialize termination flag
    termination = 0;


    % stagnation counter
    if abs(f - f_old) < opt.f_upd_tol
      i_stagnate = i_stagnate + 1;
    else
      i_stagnate = 0;
    end

    if i_stagnate >= opt.n_stagnate
    % progress has stagnated
      termination = 1;
      success = 1;      
      message = ['Optimization has stagnated (',num2str(i_stagnate),' updates with |f-f_old| < ',num2str(opt.f_upd_tol,'%.2g'),' in series)'];
    end

    if f <= opt.f_tol
    % absolute tolerance for f reached
      termination = 1;
      success = 1;      
      message = 'Absolute tolerance for f reached';
    end

    if dfnorm <= opt.dfnorm_tol
    % absolute tolerance for dfnorm reached      
      termination = 1;
      success = 1;
      message = 'Absolute tolerance for ||df|| reached.';
    end

    if dx_update <= opt.dx_upd_tol
    % too little progress in solution
      termination = 1;
      success     = 1;
      message     = 'Absolute tolerance for ||dx_update|| reached.';
    end

    if df_update <= opt.df_upd_tol
    % too little progress in gradient
      termination = 1;
      success     = 1;
      message     = 'Absolute tolerance for ||df_update|| reached.';
    end

    if iter >= opt.maxiter
    % maximum number of steps reached      
      termination = 1;
      success = 0;
      message = 'Maximum number of steps reached.';
    end

  % determine continue_itertion flag
    if termination == 1
      continue_itertion = 0;
    else
      continue_itertion = 1;
      success = NaN; % dummy value
    end

  % print info when loop terminates
    if continue_itertion == 0
      if opt.print_info > 0
        fprintf(1,'Optimization has stopped after %d iterations: ', iter)
        fprintf(1,'%s\n',message);
      end
    end
    

end