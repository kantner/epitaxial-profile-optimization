function C4 = coeff_C4(n_range, C2)
% compute bandstructure coefficients C4 via
%   C4(k) = sum_n C2(n-k)^* C2(n)

% allocate memory
  C4 = zeros(size(n_range));

% loop over k 
  for ik = 1 : length(n_range)

  % get value of k  
    k = n_range(ik);

  % loop over n
    for in = 1 : length(n_range)
    % get value of n
      n = n_range(in);
      
    % find idx in list corresponding to n-k  
      idx_1 = find(n_range == n-k);
      if idx_1>0 % proceed only if n-k is in n_range
        
        %idx_2 = find(n_range == n); % thus should always give idx_2 = in
        idx_2 = in;
        %if idx_2>0
          C4(ik) = C4(ik) + C2(idx_1)' * C2(idx_2);
        %end
      end
     
    end  
  end


end