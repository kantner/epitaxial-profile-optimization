function [] = optimization_map_scan_eval_1D(results, par)

  % determine sweep direction
    nx = length(results.x_mean_range)
    nk = length(results.kc_range)

    if nk == 1
      x_range = results.x_mean_range;
      x_title = 'mean Ge concentration x';
      fixed_val_info = ['k_c = ',num2str(results.kc_range * par.a0/(2*pi),'%.3f'), ' (2\pi/a_0)'];

    elseif nx == 1
      x_range = results.kc_range * par.a0/(2*pi);
      x_title = 'filter cutoff k_c (2\pi/a_0)';
      fixed_val_info = ['x = ',num2str(results.x_mean_range ,'%.3f')];      
    end
  


  figure(66666); clf; hold all;

  sgtitle(['linescan at ',fixed_val_info])

  %%%%%%%%%%%%%%%%%%%%%%%%  
  % plot profile
  subplot(1,3,1); hold all;
    for i = 1 : length(x_range)
  
      tmp = results.map(i).x_opt;
      tmp = apply_filter(tmp, results.map(i).par.opt.filter);
      tmp = apply_window(tmp, results.map(i).par.opt.window);
  
      tmp = par.X_QW + tmp;
  
      plot(par.z/par.units.nm, i + par.X_QW/par.X_barrier,'-','Color',[1 1 1]*0.5,'LineWidth',1)
      plot(par.z/par.units.nm, i + tmp/par.X_barrier,'k-','LineWidth',2)
      
    end
    xlabel('z (nm)')
    ylabel('profile X')
    title('epitaxital profile')
    box on;

    %%%%%%%%%%%%%%%%
    % plot mean and std
  
    subplot(1,3,2); hold all;
    M = zeros(length(x_range),1);
    S = zeros(length(x_range),1);
    for i = 1 : length(x_range)
      M(i) = results.map(i).storage.mean_EVS_list(end);
      S(i) = results.map(i).storage.std_EVS_list(end);
    end
    plot(x_range,M/par.units.meV,'ko-','DisplayName','mean E_{VS}')
    plot(x_range,S/par.units.meV,'ro-','DisplayName','std E_{VS}')
    xlabel(x_title)
    ylabel('energy (meV)')
    title('mean(E_{VS}) and std(E_{VS})')
    legend('Location','best')
    box on;
  
    subplot(1,3,3); hold all;
    plot(x_range,S./M,'ko-','DisplayName','std E_{VS}/mean E_{VS}')
    xlabel(x_title)
    ylabel('S/M')
    title('ratio std/mean')
    box on;

end