function [idx] = find_minusG(G_list)
  N_G = size(G_list,2);
  idx = zeros(1,N_G);
  for i = 1 : N_G
    idx(i) = find( sum(abs(G_list+G_list(:,i))) == 0);
  end

end