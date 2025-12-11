function [x_opt, stats, storage] = optimization(x, par, opt)

  switch(par.opt.method)
    case 1 % gradient descent
      [x_opt, stats, storage] = gradientDescent(x, par, opt);
    case 2 % Barzilai-Borwein
      [x_opt, stats, storage] = BarzilaiBorwein(x, par, opt);
    case 3 % L-BFGS
      [x_opt, stats, storage] = LBFGS(x, par, opt);
    otherwise
      error('not implemented')
  end

end