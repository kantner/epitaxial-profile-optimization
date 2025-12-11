function [] = script_Figs_2_5_6(par)
  % load to path
    addpath('optimization/')
    addpath('src/opt/')
    addpath('optimization/functions')
    addpath('optimization/user')

  % numerics
    maxiter_BB    = 0;
    maxiter_LBFGS = 100;
    plot_skip     = 100;
    restarts      = 2; % restarts of L-BFGS
    refinements   = 2; % ramp up w(2) weight after restarts to meet constraint
    continuation  = 0; % 0 = start from zero | 1 = from previous profile   
    conventional  = 0; % 0 = off | 1 = skip optimization, make plot with conventional WW that has same Ge budget and q=kc

  % physics
    par.eps      = par.eps_QW;
    par.eps(1,2) = 0.1 * par.units.percent;
    par.eps(2,1) = 0.1 * par.units.percent;    
    F = 5*1E6;

  % data storage for comparison    
    compare_opt = [];

 %% optimization for useful properties
  % case A --> min 1/nu
    optimization_objective = 13;
    X_budget  = 0.05;
    kc_list   = [0.5] * 2*pi/par.a0;     
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);

    compare_opt(1).name  = 'modulated wiggle well (A)';
    compare_opt(1).out   = out;
    compare_opt(1).x_opt = x_opt;
    compare_opt(1).x_mod = x_mod;

  % case B --> min sqrt(2*Gamma)/nu
    optimization_objective = 12;
    X_budget  = 0.05;
    kc_list   = [0.5] * 2*pi/par.a0;           
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);
    
    compare_opt(2).name  = 'modulated wiggle well (B)';
    compare_opt(2).out   = out;
    compare_opt(2).x_opt = x_opt;
    compare_opt(2).x_mod = x_mod;

  % case C --> min sqrt(2*Gamma)
    optimization_objective = 14;
    X_budget  = 0.05;
    kc_list   = [0.5] * 2*pi/par.a0;           
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);
 
    compare_opt(3).name  = 'narrow well';
    compare_opt(3).out   = out;
    compare_opt(3).x_opt = x_opt;
    compare_opt(3).x_mod = x_mod;

 %  Ge spike
    optimization_objective = 13; % min 1/nu
    X_budget  = 0.05;
    kc_list   = 0.07 * 2*pi/par.a0;     
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par)            

    compare_opt(4).name  = 'Ge spike';
    compare_opt(4).out   = out;
    compare_opt(4).x_opt = x_opt;
    compare_opt(4).x_mod = x_mod;

 %% profiles for reference
  % conventional WW
    kc_list      = 2*par.k1;
    X_budget     = 0.05;
    conventional = 1;
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);               
    conventional = 0;  % switch back

    compare_opt(5).name  = 'conventional wiggle well';
    compare_opt(5).out   = out;
    compare_opt(5).x_opt = x_opt;
    compare_opt(5).x_mod = x_mod;

  % plain smoothed QW
    kc_list      = 1.0;
    X_budget     = 0.0;
    conventional = 1;
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);               
    conventional = 0;  % switch back

    compare_opt(6).name  = 'plain (smoothed) QW';
    compare_opt(6).out   = out;
    compare_opt(6).x_opt = x_opt;
    compare_opt(6).x_mod = x_mod;

  % plain sharp QW
    kc_list      = 2.0;
    X_budget     = 0.0;

    par.sigma_l = 1E-5 * par.units.nm;
    par.sigma_u = 1E-5 * par.units.nm;
    par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+0.5*par.h_QW)/par.sigma_u);  % QW at [-h/2 h/2]
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);
    par.U_QW = par.dEc * par.X_QW;

    conventional = 1;
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);               
    conventional = 0;  % switch back

    compare_opt(7).name  = 'plain (sharp) QW';
    compare_opt(7).out   = out;
    compare_opt(7).x_opt = x_opt;
    compare_opt(7).x_mod = x_mod;

    % switch back
    par.sigma_l = 0.5 * par.units.nm;
    par.sigma_u = 0.5 * par.units.nm;
    par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+0.5*par.h_QW)/par.sigma_u);  % QW at [-h/2 h/2]
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);
    par.U_QW = par.dEc * par.X_QW;


 %% optimization for inverted objectives
  % case D (opposite of A) --> minimize nu (flat Ge)
    optimization_objective = 15;
    X_budget = 0.05;
    kc_list       = [0.5] * 2*pi/par.a0;     
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par)

    compare_opt(8).name  = 'flat Ge';
    compare_opt(8).out   = out;
    compare_opt(8).x_opt = x_opt;
    compare_opt(8).x_mod = x_mod;    

  % case E (opposite of C) --> maximize Gamma (Ge plateau)
    optimization_objective = 16; 
    X_budget = 0.05;
    kc_list       = [0.5] * 2*pi/par.a0;     
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par)      

    compare_opt(9).name  = 'max disorder';
    compare_opt(9).out   = out;    
    compare_opt(9).x_opt = x_opt;
    compare_opt(9).x_mod = x_mod;

 %% case F --> maximize variance
    optimization_objective = 4; 
    X_budget = 0.05;
    kc_list       = [0.5] * 2*pi/par.a0;     
    [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par)      

    compare_opt(10).name  = 'max variance';
    compare_opt(10).out   = out;    
    compare_opt(10).x_opt = x_opt;
    compare_opt(10).x_mod = x_mod;    
        


 %% make figure for comparison
    marker_symbol_list = {'d','*','o','v','s','d','+','o','v','s'}
    cmap = lines(12);

    fig_comparison = figure(13425);clf; hold all;      
    for i = 1 : length(compare_opt)
      plot(compare_opt(i).out.nu/par.units.ueV, sqrt(2*compare_opt(i).out.Gamma/par.units.ueV^2), 'o', 'DisplayName',compare_opt(i).name, 'MarkerFaceColor',cmap(i,:),'Marker',marker_symbol_list{i},'Color',cmap(i,:),'MarkerSize',10)
    end
    legend('Location','best')
    box on;
    xlabel('deterministic \nu (ueV)')
    ylabel('random \sqrt(2 \Gamma) (ueV)')

    nu_range = linspace(0,1000,101);
    plot(nu_range, 0.3507 * 2*nu_range, 'k--')

    xlim([-25 1025])
    ylim([-5 105])

    exportgraphics(fig_comparison,'comparison_map.pdf','ContentType','vector')


end