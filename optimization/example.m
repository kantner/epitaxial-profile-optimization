clear all
clc
profile on
%% add folders
   addpath('functions/')
   addpath('user/')

%% Example for optimization
   % random number seed
     rng(12394)
    
  % parameter
    par.model = 1; % 1 = Rosenbrock | 2 = Styblinski–Tang
    switch(par.model)
      case 1
        par.a = 1;
        par.b = 100;
      case 2
        par.a = -16;
        par.b = 5;
    end
    par.N = 2;
    
  % initialization  
    x0 = randn(par.N,1);
    %x0 = [-1.5; 2];
    
  % prepare map
    opt.map.x     = linspace(-2,2,101);
    opt.map.y     = linspace(-1,3,102);

    opt.map.x     = linspace(-5,5,101);
    opt.map.y     = linspace(-5,5,102);
    opt.map.f_map = f_map_2D(par, opt.map.x, opt.map.y);
    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % optimization method
    opt.method     = 1;     % 1 = gradient descend | 2 = Barzilai-Borwein | 3 = L-BFGS
    opt.maxiter    = 10000; 
    opt.f_tol      = -1E13;
    opt.dfnorm_tol = 1E-13;
    opt.print_info = 1; % 0 = off | 1 = basic info | 2 = detailed info
    opt.plot_skip  = 1000; % update plot after plot_skip steps 
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % method specific parameters

  % Armijo linesearch parameters
    opt.Armijo.alpha0  = 1;
    opt.Armijo.maxiter = 25;
    opt.Armijo.c1      = 0.5;
    opt.Armijo.shrink  = 0.5;
    opt.Armijo.increase = 1.5;
    opt.Armijo.print_info = 0; % 0 = off | 1 = on
    
  % Wolfe linesearch parameters
    opt.Wolfe.alpha0  = 1;
    opt.Wolfe.maxiter = 25;
    opt.Wolfe.c1      = 0.5;
    opt.Wolfe.c2      = 0.999;
    opt.Wolfe.lambda  = 2;
    opt.Wolfe.print_info = 0; % 0 = off | 1 = on
  
  % gradient descent parameters
    opt.gradientDescent.linesearch_mode = 1; % 1 = Armijo | 2 = Wolfe
    
  % Barzilai Borwein parameters
    opt.BB.alpha0     = 1;
    opt.BB.alpha_mode = 6; % mode for step size computation: 1 = long | 2 = short | 3 = geom. mean | 4 = alternating | 5 = random | 6 = positive alternating
    
  % L-BFGS parameters
    opt.LBFGS.memory = 20;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % run optimization  
    [x] = gradientDescent(x0, par, opt);
    pause
    [x] = BarzilaiBorwein(x0, par, opt);
pause
    [x] = LBFGS(x0, par, opt);

    profile viewer


   