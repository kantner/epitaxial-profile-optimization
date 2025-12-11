function [a] = scalar_product(p, df, par)
  a = par.dz * (p'*df);
end
