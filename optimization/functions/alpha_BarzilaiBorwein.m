function alpha = alpha_BarzilaiBorwein(s,y,iter,opt)
% stepsize selection for Barzilai-Borwein method

% check for errors
  if sum(isnan(y))>0
    error('y contains NaN')
  end

  if sum(isinf(y))>0
    error('y contains Inf')
  end

  if sum(isnan(s))>0
    error('s contains NaN')
  end

  if sum(isinf(s))>0
    error('s contains Inf')
  end

  if sum(abs(y)) == 0
    error('y identical zero')
  end

% alpha long and short  
  s2 = s'*s;
  y2 = y'*y;
  sy = s'*y;

  alpha_long  = s2/sy;
  alpha_short = sy/y2;

% set step size  
  switch(opt.BB.alpha_mode)
    case 1 % long
      alpha = alpha_long;
    case 2 % short
      alpha = alpha_short;
    case 3 % geometric mean
      alpha = sqrt(alpha_long * alpha_short);
    case 4 % alternating
      if mod(iter,2) == 0
        alpha = alpha_long;        
      else
        alpha = alpha_short;
      end
    case 5 % random
      u = rand;
      if u > 0.5
        alpha = alpha_long;        
      else
        alpha = alpha_short;
      end
    case 6 % positive alternating
      % if both alpha are positive: alternating
      % else: geometric mean
      if sy > 0
        if mod(iter,2) == 0
          alpha = alpha_long;        
        else
          alpha = alpha_short;
        end
      else
        alpha = sqrt(alpha_long * alpha_short);
      end


    otherwise
      error('not implemented')
  end


% checks
  if isnan(alpha)
    error('alpha is NaN')
  end

  if isinf(alpha)
    error('alpha is Inf')
  end


end