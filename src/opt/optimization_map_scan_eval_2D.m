function [] = optimization_map_scan_eval_2D(results, par)


  figure(66666); clf; hold all;

%% EVAL 2D MAP SCAN
  subplot(2,2,1); hold all;
  for ik = 1 : length(results.kc_range)
    for im = 1 : length(results.x_mean_range)

    tmp = results.map(ik,im).x_opt;
    tmp = apply_filter(tmp, results.map(ik,im).par.opt.filter);
    tmp = apply_window(tmp, results.map(ik,im).par.opt.window);

    tmp = par.X_QW + tmp;

    plot(par.z/par.units.nm, ik + tmp/par.X_barrier,'k-')
    end
  end
  xlabel('z (nm)')
  ylabel('profile X')
  title('epitaxital profile')
  
  
  M = zeros(length(results.kc_range),length(results.x_mean_range));
  S = zeros(length(results.kc_range),length(results.x_mean_range));
  for ik = 1 : length(results.kc_range)
    for im = 1 : length(results.x_mean_range)
    M(ik,im) = results.map(ik,im).storage.mean_EVS_list(end);
    S(ik,im) = results.map(ik,im).storage.std_EVS_list(end);
    end
  end

  cmap = divergingColormap1(1024);
  subplot(2,2,2); hold all;
    surf(results.kc_range * par.a0/(2*pi), results.x_mean_range, (M)'/par.units.meV)
    xlabel('k_c (2\pi/a_0)')
    ylabel('x mean')
    title('mean')
    box on
    colorbar
    %shading interp
    colormap(cmap)

  subplot(2,2,3); hold all;
    surf(results.kc_range * par.a0/(2*pi), results.x_mean_range, (S)'/par.units.meV)
    xlabel('k_c (2\pi/a_0)')
    ylabel('x mean')
    title('std')
    box on
    colorbar
    colormap(cmap)


  subplot(2,2,4); hold all;
    surf(results.kc_range * par.a0/(2*pi), results.x_mean_range, (S./M)')
    xlabel('k_c (2\pi/a_0)')
    ylabel('x mean')
    title('std/mean')
    box on
    colorbar
    colormap(cmap)
    %shading interp
    %set(gca,'ZScale','log')
    %set(gca,'ColorScale','log')

end