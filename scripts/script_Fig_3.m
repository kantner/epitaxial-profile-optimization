function [] = script_Fig_3(par)
    
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


  % run sweep  
    X_budget_list = [0.01 : 0.01 : 0.15];
    for iX = 1 : length(X_budget_list)
      optimization_objective = 12; % --> min sqrt(2*Gamma)/nu
      kc_list       = [0.5] * 2*pi/par.a0;     
      X_budget = X_budget_list(iX);
      [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par)
      mod_WW_B(iX).out   = out;
      mod_WW_B(iX).x_opt = x_opt;
      mod_WW_B(iX).x_mod = x_mod;
    end

    %%
  fig_Ge_sweep = figure(3);clf; hold all;
    
  subplot(1,2,1); hold all; 
  cmap = dusk(13);
  
  plot(par.z/par.units.nm, par.X_QW,'k--','DisplayName','0%','Color','k')
  %i = 1;
  %plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'k-','DisplayName','1%','Color',cmap(i,:))
  i = 2;
  plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'k-','DisplayName','2%','Color',cmap(i,:))
  %i = 3;
  %plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'r-','DisplayName','3%','Color',cmap(i,:))
  i = 4;
  plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','4%','Color',cmap(i,:))
  %i = 5;
  %plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','5%','Color',cmap(i,:))
  i = 6;
  plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','6%','Color',cmap(i,:))
  %i = 7;
  %plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','7%','Color',cmap(i,:))
  i = 8;
  plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','8%','Color',cmap(i,:))
  %i = 9;
  %plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','9%','Color',cmap(i,:))
  i = 10;
  plot(par.z/par.units.nm, par.X_QW + mod_WW_B(i).x_mod,'g-','DisplayName','10%','Color',cmap(i,:))
  box on
  legend
  xlim([-11,11])
  ylim([-0.025,0.375])
  xlabel('space z (nm)')
  ylabel('Ge conc. ')

  subplot(1,2,2); hold all; 
  for i = 1 : 13%;length(mod_WW_B)
    plot(mod_WW_B(i).out.nu/par.units.ueV,sqrt(2*mod_WW_B(i).out.Gamma)/par.units.ueV,'ko','MarkerSize',10,'Color',cmap(i,:),'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','k')
  end
  xlim([-50,2050])
  ylim([-5,105])
  xlabel('\nu (ueV)')
  ylabel('\sqrt{2*\Gamma} (ueV)')
  box on


  exportgraphics(fig_Ge_sweep,'mod_WW_B_Ge_sweep.pdf','ContentType','vector')


end