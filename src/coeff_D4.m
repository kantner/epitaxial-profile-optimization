function D4 = coeff_D4(n_range, C2)
% compute bandstructure coefficients D4 via
%   D4(n) = sum_m C2(m) C2(n-m)

% allocate memory
  D4 = zeros(size(n_range));

% loop over m 
  for im = 1 : length(n_range)

  % get value of m  
    m = n_range(im);

  % loop over n
    for in = 1 : length(n_range)
    % get value of n
      n = n_range(in);
      
    % find idx in list corresponding to n-m
      idx_2 = find(n_range == n-m);
      if idx_2>0 % proceed only if n-k is in n_range
        
        idx_1 = im; 
        
        D4(in) = D4(in) + C2(idx_1) * C2(idx_2);

      end
     
    end  
  end


end