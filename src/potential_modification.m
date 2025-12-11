function U_x = potential_modification(x, par)

  %U_x = par.dEc * par.opt.window .* x;
  U_x = par.dEc .* x;

end