function [f,df] = f_penalty(x, par)
% penality function to enforce x in admissible range a<=x<=b

% import
  eps = par.opt.penalty.eps;
  b   = par.opt.penalty.x_max;
  a   = par.opt.penalty.x_min;

% check for errors  
  if eps<0
    error('eps must be non-negative!')
  end

  if b<a
    error('b must be larger than a')
  end

% evaluate  
  if eps > 0
  % function
    % f = 0.5*eps*log(4*cosh((x-a)/eps).*cosh((x-b)/eps)) - 0.5*abs(b-a);
    % this above expression is not robust for numerical evaluation
  
  % robust evaluation  
    f =  - 0.5*abs(b-a) ...
         + 0.5*eps*(log(4) +  log_cosh((x-a)/eps) + log_cosh((x-b)/eps) );
     
  % derivative 
    if nargout == 2
      df = 0.5 * ( tanh((x-a)/eps) + tanh((x-b)/eps) );
    end


  else % eps = 0

    % piecewise definition
      f = zeros(size(x));
    
      idx = x>b;
      f(idx) = x(idx) - b;
    
      idx = x < a;
      f(idx) = a - x(idx);
      

    % derivative 
    if nargout == 2
      
      df = zeros(size(x));
    
      df(x > b) = +1;    
      df(x < a) = -1;

      
    end
  
  end

end