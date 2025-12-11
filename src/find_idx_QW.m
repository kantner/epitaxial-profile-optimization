function [idx] = find_idx_QW(par)

  g_list = par.G_list*par.a0/(2*pi);
  

  tol = 1E-12;
  idx = [];
  for i = 1 : par.N_G
    gi = g_list(:,i);
    for j = 1 : par.N_G
      gj = g_list(:,j);
      if abs(gi(1)+gj(1)) + abs(gi(2)+gj(2)) <= tol
        % add if pair satisfies selection rule
        idx = [idx; i, j];
      end
    end

    

  end


end