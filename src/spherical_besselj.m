function [j] = spherical_besselj(n,x)
  
%{
  if x == 0

    if n == 0
      j = 1;
    else
      j = 0;
    end
  else
    j = sqrt(0.5*pi/x) .* besselj(n+0.5,x);
  end


%}


switch(n)
  case -1
     j = cos(x)./x;    

  case 0
    if x == 0
      j = 1;%1 - (1/6) * x*x;% + (1/120) * x(idx).^4;
    else
      j = sin(x)./x;
    end
    
  case 1
    if x == 0
      j = 0;%1/3 * x .*( 1 - 1/10 * x.*x);
    else
      j = (sin(x)./x - cos(x))./x;
    end

  otherwise
    error('not implemented')

end


end