function [] = print_opt_progress(iter, f, dfnorm, f_old, opt)
% plot optimization progress

  if opt.print_info >= 2
    
    % print headline
    if mod(iter,20) == 1
      fprintf(1,'iter\tf(x)\t\t||df(x)||\t|f-f_old|\n')
    end

    % print current values
      fprintf(1,'%d\t%.4e\t%.4e\t%.4e\n',iter, f, dfnorm, abs(f-f_old));

  end

end