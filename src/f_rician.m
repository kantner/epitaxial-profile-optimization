function [f, df] = f_rician(x)

  I0 = besseli(0, x);
  I1 = besseli(1, x);
  expmx = exp(-x);

%%%%%%%%%%%%%
% function f
  f = expmx.*( (1+2*x).*I0 + 2*x .* I1);

% asymptotics for large values
  x_high = 700;
  idx = x>x_high;

  f(idx) = 2*sqrt(2*x(idx)/pi) + 1./(2*sqrt(2*pi*x(idx)))+ 1./(32.*sqrt(2*pi).*x(idx).^1.5 );


  
%%%%%%%%%%%%%
% derivative  
  if nargout == 2
    df      = expmx.*(I0 + I1);

    % asymptotics for large values  
    df(idx) = sqrt(2/pi * 1./x(idx)) - 1./( 4 * sqrt(2*pi) * x(idx).^1.5);
  end

end