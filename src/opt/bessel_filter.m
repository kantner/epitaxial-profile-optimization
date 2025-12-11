function [F] = bessel_filter(n, omega, omega_c)

  % coefficients
    a = zeros(n+1);

    for k = 0 : n
      a(k+1) = factorial(2*n-k)/(2^(n-k) * factorial(k) * factorial(n-k));
    end


  % numerator
    numer = a(1);

  % build ap denominator (unscaled)
    denom = @(s) a(1);
    for k = 1 : n
      denom = @(s) denom(s) + a(k+1)*s.^k;
    end

  % find root to scale frequency
    scale_factor = fzero(@(s) denom(s)-2*numer,0.7);

  % generate scaled Bessel filter  
    F = numer./denom( scale_factor * abs(omega)/omega_c  );


    

end