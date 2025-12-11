function [results] = optimization_map_scan(x_mean_range, kc_range, par)

  results.x_mean_range  = x_mean_range;
  results.kc_range      = kc_range;

  % initialize
    x = zeros(par.N,1);    

  for ik = 1 : length(kc_range)
    for im = 1 : length(x_mean_range)

    ik
    im


  % initialize
    x = zeros(par.N,1);            

  % set parameters
    par.opt.k_c      = kc_range(ik);
    par.opt.X_budget = x_mean_range(im);

  % update filter
    par.opt.filter       = generate_filter(par.k, par.opt.k_c, par.opt.filter_type);

  % run optimization: stage 1 (BB)
  %
    par.opt.plot_skip = 100; % update plot after plot_skip steps   
    par.opt.maxiter = 500; 
    par.opt.method  = 2; 
    [x_opt, stats, storage] = optimization(x, par, par.opt);

  % transfer initialization  
    x = x_opt;
  %{
    x = apply_filter(x, par.opt.filter);
    x = apply_window(x, par.opt.window);
  %}

  % run optimization: stage 2 (L-BFGS)   
  %
    par.opt.LBFGS.memory = 1000;
    par.opt.plot_skip = 100; % update plot after plot_skip steps   
    par.opt.maxiter   = 1000; 
    par.opt.method    = 3; 
    [x_opt, stats, storage] = optimization(x, par, par.opt);

  % filter and window results  
  %
    x = x_opt;
    x = apply_filter(x, par.opt.filter);
    x = apply_window(x, par.opt.window);
  %}
  % store
    results.map(ik,im).x_opt   = x_opt;
    results.map(ik,im).storage = storage;
    results.map(ik,im).par     = par;

  % update initialization
    x = x;
  

    disp('---- done ----')  
    %pause

    end


  end