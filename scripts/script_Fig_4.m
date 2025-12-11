function [] = script_Fig_4(par)
  % numerics
    maxiter_BB    = 0;
    maxiter_LBFGS = 100;
    plot_skip     = 100;
    restarts      = 2; % restarts of L-BFGS
    refinements   = 2; % ramp up w(2) weight after restarts to meet constraint
    continuation  = 0; % 0 = start from zero | 1 = from previous profile   
    conventional  = 0; % 0 = off | 1 = skip optimization, make plot with conventional WW that has same Ge budget and q=kc

  % physics
    physical_constants;
    par.eps      = par.eps_QW;
    par.eps(1,2) = 0.1 * par.units.percent;
    par.eps(2,1) = 0.1 * par.units.percent;    

  % design field strength  
    F = 5*1E6;

  % compute modulated wiggle well case B --> mod WW (for given field F)
    optimization_objective = 12;
    kc_list       = [0.4] * 2*pi/par.a0;     
    X_budget = 0.05;

    [x_opt, x_modWW, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);

  % compute Ge spike
    optimization_objective = 13;
    kc_list       = [0.07] * 2*pi/par.a0;     
    X_budget = 0.05;

    [x_opt, x_Gespike, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par);
  
  % conventional WW for comparison
    A = compute_wiggle_well_amplitude(X_budget, 2*par.k1,0,par);
    x_convWW = x_wiggle_well(A, 2*par.k1, 0, par);


 %% sweep over different field strength
    F_list = unique([-5,0,5,10,linspace(-13,13,350)]) * 1E6; % mV/nm

    for iF = 1 : length(F_list)
      par.F    = F_list(iF);
      par.U_F  = -elementaryCharge*par.F * par.z;
      res_modWW(iF).out   = compute_valley_splitting(x_modWW, 0, par);
      res_convWW(iF).out  = compute_valley_splitting(x_convWW, 0, par);
      res_Gespike(iF).out = compute_valley_splitting(x_Gespike, 0, par);
    end

  % extract data for plot    
    mean_modWW  = zeros(1,length(F_list));
    mean_convWW = zeros(1,length(F_list));
    mean_Gespike = zeros(1,length(F_list));

    nu_modWW   = zeros(1,length(F_list));
    nu_convWW  = zeros(1,length(F_list));
    nu_Gespike = zeros(1,length(F_list));

    Gamma_modWW   = zeros(1,length(F_list));
    Gamma_convWW  = zeros(1,length(F_list));
    Gamma_Gespike = zeros(1,length(F_list));

    for iF = 1 : length(F_list)
      mean_modWW(iF)  = res_modWW(iF).out.M;
      mean_convWW(iF) = res_convWW(iF).out.M;
      mean_Gespike(iF) = res_Gespike(iF).out.M;

      nu_modWW(iF)   = res_modWW(iF).out.nu;
      nu_convWW(iF)  = res_convWW(iF).out.nu;
      nu_Gespike(iF) = res_Gespike(iF).out.nu;

      Gamma_modWW(iF)  = res_modWW(iF).out.Gamma;
      Gamma_convWW(iF) = res_convWW(iF).out.Gamma;
      Gamma_Gespike(iF) = res_Gespike(iF).out.Gamma;
    end

  % compute confidence intervals
    upper_bound = 0.75;
    lower_bound = 0.25;   

    sigma_lower_modWW = zeros(1,length(F_list));
    sigma_upper_modWW = zeros(1,length(F_list));

    sigma_lower_convWW = zeros(1,length(F_list));
    sigma_upper_convWW = zeros(1,length(F_list));

    sigma_lower_Gespike = zeros(1,length(F_list));
    sigma_upper_Gespike = zeros(1,length(F_list));

    for iF = 1 : length(F_list)    

      nu    = nu_modWW(iF);
      sigma = sqrt(2*Gamma_modWW(iF));
      sigma_lower_modWW(iF) = par.units.ueV * compute_rice_percentile(lower_bound, nu/par.units.ueV, sigma/par.units.ueV);
      sigma_upper_modWW(iF) = par.units.ueV * compute_rice_percentile(upper_bound, nu/par.units.ueV, sigma/par.units.ueV);

      nu    = nu_convWW(iF);
      sigma = sqrt(2*Gamma_convWW(iF));
      sigma_lower_convWW(iF) = par.units.ueV * compute_rice_percentile(lower_bound, nu/par.units.ueV, sigma/par.units.ueV);
      sigma_upper_convWW(iF) = par.units.ueV * compute_rice_percentile(upper_bound, nu/par.units.ueV, sigma/par.units.ueV);

      nu    = nu_Gespike(iF);
      sigma = sqrt(2*Gamma_Gespike(iF));
      sigma_lower_Gespike(iF) = par.units.ueV * compute_rice_percentile(lower_bound, nu/par.units.ueV, sigma/par.units.ueV);
      sigma_upper_Gespike(iF) = par.units.ueV * compute_rice_percentile(upper_bound, nu/par.units.ueV, sigma/par.units.ueV);

    end

 %% plot
    fig_F_sweep = figure(21324);clf;hold all;

    nrows = 3;
    ncols = 3; 
   
    subplot(nrows, ncols, 1); hold all;
      plot(F_list * 1E-6, mean_modWW/par.units.ueV, 'r-','DisplayName','mod WW')
      plot(F_list * 1E-6, mean_convWW/par.units.ueV, 'b-','DisplayName','conv WW')
      plot(F_list * 1E-6, mean_Gespike/par.units.ueV, 'g-','DisplayName','Ge spike')
      
      plot(F_list * 1E-6, (sigma_lower_modWW)/par.units.ueV, 'r--','DisplayName','upper conf. mod WW')
      plot(F_list * 1E-6, (sigma_upper_modWW)/par.units.ueV, 'r--','DisplayName','lower conf. mod WW')

      plot(F_list * 1E-6, (sigma_lower_convWW)/par.units.ueV, 'b--','DisplayName','upper conf. conv WW')
      plot(F_list * 1E-6, (sigma_upper_convWW)/par.units.ueV, 'b--','DisplayName','lower conf. conv WW')    
  
      plot(F_list * 1E-6, (sigma_lower_Gespike)/par.units.ueV, 'g--','DisplayName','upper conf. Ge spike')
      plot(F_list * 1E-6, (sigma_upper_Gespike)/par.units.ueV, 'g--','DisplayName','lower conf. Ge spike')
  
      xline(5,'k--','DisplayName','F_{opt}')

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('mean valley splitting \langleE_{VS}\rangle (ueV)')
      xlim([-1 1]*12)
      ylim([-75 1275])
      set(gca,'YTick',[-0:100:1200])

    subplot(nrows, ncols, 2); hold all;      
      plot(0.5*(F_list(1:end-1)+F_list(2:end)) * 1E-6, (mean_modWW(2:end)-mean_modWW(1:end-1))./(F_list(2:end)-F_list(1:end-1)) * 1E6/par.units.ueV, 'r-','DisplayName','mod WW')
      plot(0.5*(F_list(1:end-1)+F_list(2:end)) * 1E-6, (mean_convWW(2:end)-mean_convWW(1:end-1))./(F_list(2:end)-F_list(1:end-1)) * 1E6/par.units.ueV, 'b-','DisplayName','conv WW')
      plot(0.5*(F_list(1:end-1)+F_list(2:end)) * 1E-6, (mean_Gespike(2:end)-mean_Gespike(1:end-1))./(F_list(2:end)-F_list(1:end-1)) * 1E6/par.units.ueV, 'g-','DisplayName','Ge spike')

      xline(5,'k--','DisplayName','F_{opt}')
      box on
      legend('Location','best')      
      xlabel('electric field F (mV/nm)')
      ylabel('field tunability \partial\langleE_{VS}\rangle/\partial F (ueV*nm/mV)')
      xlim([-1 1]*12)
      ylim([-70 470])
      set(gca,'YTick',[-50:50:450])


    subplot(nrows, ncols, 3); hold all;
      plot(F_list * 1E-6, 0.5*nu_modWW./sqrt(2*Gamma_modWW), 'r-','DisplayName','mod WW')
      plot(F_list * 1E-6, 0.5*nu_convWW./sqrt(2*Gamma_convWW), 'b-','DisplayName','conv WW')
      plot(F_list * 1E-6, 0.5*nu_Gespike./sqrt(2*Gamma_Gespike), 'g-','DisplayName','conv WW')
  
      
      xline(F*1E-6,'k--','DisplayName','F_{opt}')

      box on
      legend('Location','best')         
      xlabel('electric field F (mV/nm)')
      ylabel('deterministic fraction factor \zeta')
      yline(0.3507,'m:','DisplayName','Q=1/2')
      xlim([-1 1]*12)      
      ylim([-0.5 6.5])
      

    subplot(nrows, ncols, 4); hold all;
      F_target = -5 * 1E6;
      idx = find(abs(F_list - F_target) == min(abs(F_list - F_target)));

      plot(par.z/par.units.nm, (-elementaryCharge*par.z*F_list(idx) + par.dEc*(x_modWW+par.X_QW))/par.units.meV,'k-','DisplayName',['F = ',num2str(F_list(idx)*1E-6),' mV/nm'])
      plot(par.z/par.units.nm, res_modWW(idx).out.E(res_modWW(idx).out.idx_gnd)/par.units.meV + par.N * (res_modWW(idx).out.S(:,res_modWW(idx).out.idx_gnd)).^2,'r-')         

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('energy (meV)')
      maxU = 265;
      minU = -65;
      xlim([-11,11])
      ylim([minU,maxU])


    subplot(nrows, ncols, 5); hold all;
      F_target = +5 * 1E6;
      idx = find(abs(F_list - F_target) == min(abs(F_list - F_target)));
      plot(par.z/par.units.nm, (-elementaryCharge*par.z*F_list(idx) + par.dEc*(x_modWW +par.X_QW))/par.units.meV,'k-','DisplayName',['F = ',num2str(F_list(idx)*1E-6),' mV/nm'])
      plot(par.z/par.units.nm, res_modWW(idx).out.E(res_modWW(idx).out.idx_gnd)/par.units.meV + par.N * (res_modWW(idx).out.S(:,res_modWW(idx).out.idx_gnd)).^2,'r-')

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('energy (meV)')
      xlim([-11,11])
      ylim([minU,maxU])

  subplot(nrows, ncols, 6); hold all;
      F_target = +10 * 1E6;
      idx = find(abs(F_list - F_target) == min(abs(F_list - F_target)));
      plot(par.z/par.units.nm, (-elementaryCharge*par.z*F_list(idx) + par.dEc*(x_modWW+par.X_QW))/par.units.meV,'k-','DisplayName',['F = ',num2str(F_list(idx)*1E-6),' mV/nm'])
      plot(par.z/par.units.nm, res_modWW(idx).out.E(res_modWW(idx).out.idx_gnd)/par.units.meV + par.N * (res_modWW(idx).out.S(:,res_modWW(idx).out.idx_gnd)).^2,'r-')

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('energy (meV)')
      xlim([-11,11])
      ylim([minU,maxU])




    subplot(nrows, ncols, 7); hold all;
      F_target = -5 * 1E6;
      idx = find(abs(F_list - F_target) == min(abs(F_list - F_target)));

      plot(par.z/par.units.nm, (-elementaryCharge*par.z*F_list(idx) + par.dEc*(x_Gespike+par.X_QW))/par.units.meV,'k-','DisplayName',['F = ',num2str(F_list(idx)*1E-6),' mV/nm'])
      plot(par.z/par.units.nm, res_Gespike(idx).out.E(res_Gespike(idx).out.idx_gnd)/par.units.meV + par.N * (res_Gespike(idx).out.S(:,res_Gespike(idx).out.idx_gnd)).^2,'r-')         

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('energy (meV)')
      maxU = 265;
      minU = -65;
      xlim([-11,11])
      ylim([minU,maxU])


    subplot(nrows, ncols, 8); hold all;
      F_target = +5 * 1E6;
      idx = find(abs(F_list - F_target) == min(abs(F_list - F_target)));
      plot(par.z/par.units.nm, (-elementaryCharge*par.z*F_list(idx) + par.dEc*(x_Gespike +par.X_QW))/par.units.meV,'k-','DisplayName',['F = ',num2str(F_list(idx)*1E-6),' mV/nm'])
      plot(par.z/par.units.nm, res_Gespike(idx).out.E(res_Gespike(idx).out.idx_gnd)/par.units.meV + par.N * (res_Gespike(idx).out.S(:,res_Gespike(idx).out.idx_gnd)).^2,'r-')

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('energy (meV)')
      xlim([-11,11])
      ylim([minU,maxU])

  subplot(nrows, ncols, 9); hold all;
      F_target = +10 * 1E6;
      idx = find(abs(F_list - F_target) == min(abs(F_list - F_target)));
      plot(par.z/par.units.nm, (-elementaryCharge*par.z*F_list(idx) + par.dEc*(x_Gespike+par.X_QW))/par.units.meV,'k-','DisplayName',['F = ',num2str(F_list(idx)*1E-6),' mV/nm'])
      plot(par.z/par.units.nm, res_Gespike(idx).out.E(res_Gespike(idx).out.idx_gnd)/par.units.meV + par.N * (res_Gespike(idx).out.S(:,res_Gespike(idx).out.idx_gnd)).^2,'r-')

      box on
      legend('Location','best')
      xlabel('electric field F (mV/nm)')
      ylabel('energy (meV)')
      xlim([-11,11])
      ylim([minU,maxU])


      exportgraphics(fig_F_sweep,'mod_F_sweep.pdf','ContentType','vector')

end
