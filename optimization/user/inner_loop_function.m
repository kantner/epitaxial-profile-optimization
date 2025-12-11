function [storage] = inner_loop_function(iter, alpha, x, f, df, p, out, par, opt, storage, continue_itertion)
% Function that is executed in each optimzation step.
% This function can be used to compute auxiliary quantities in a struct (storage)
% to control progress or to create intermediate plots.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize
    if iter == 1      
      storage.J_tot_list = zeros(opt.maxiter,1);
      storage.J0_list    = zeros(opt.maxiter,1);
      storage.J1_list    = zeros(opt.maxiter,1);
      storage.J2_list    = zeros(opt.maxiter,1);
      storage.J3_list    = zeros(opt.maxiter,1);
      storage.DJ_sq_list = zeros(opt.maxiter,1);
      storage.p_DJ_list  = zeros(opt.maxiter,1);
      storage.alpha_hist = zeros(opt.maxiter,1);

      storage.X_mean_1   = zeros(opt.maxiter,1);
      storage.X_mean_2   = zeros(opt.maxiter,1);

      storage.nu_list    = zeros(opt.maxiter,1);
      storage.sigma_list = zeros(opt.maxiter,1);
    end

  %%%%%%%%%%%%%%%%%%%%%
  % apply filter and window
    if isfield(par,'opt')
    % apply filter
      if par.opt.apply_filter == 1
        [x_mod] = apply_filter(x, par.opt.filter);
      end
    % apply window
      if par.opt.apply_window == 1
        [x_mod] = apply_window(x_mod, par.opt.window);
      end
    else
      x_mod = x;
    end  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % store values for progress
    storage.J_tot_list(iter)   = out.J_tot;
    storage.J0_list(iter)      = out.J0;
    storage.J1_list(iter,:)    = out.J1;
    storage.J2_list(iter,:)    = out.J2;
    storage.J3_list(iter,:)    = out.J3;
    storage.DJ_sq_list(iter,:) = dfnorm_function(out.DJ,par);
    storage.p_DJ_list(iter,:)  = -scalar_product(p, out.DJ, par);
    storage.alpha_hist(iter)   = alpha;

    %fprintf('do we need to use filtered x here?\n')
    %storage.X_mean_1(iter) = mean_Ge_budget(1, x, out.psi0, par);  
    %storage.X_mean_2(iter) = mean_Ge_budget(2, x, out.psi0, par);  
    storage.X_mean_1(iter) = mean_Ge_budget(1, x_mod, out.psi0, par);  
    storage.X_mean_2(iter) = mean_Ge_budget(2, x_mod, out.psi0, par);  

    storage.nu_list(iter)       = out.nu;
    storage.sigma_list(iter)    = out.sigma;
    storage.mean_EVS_list(iter) = out.M;
    storage.std_EVS_list(iter)  = sqrt(out.V);


  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot
    if or(mod(iter,opt.plot_skip)==0, continue_itertion == 0) % either after plot_skip or when terminating
      figure(1); clf; hold all;
        nrows   = 3;
        ncols   = 5;
        counter = 0;
    
      switch(opt.method)
        case 1
          opt_method = 'gradient descent';
        case 2
          opt_method = 'Barzilai-Borwein';
        case 3
          opt_method = 'L-BFGS';
        otherwise
          error('not implemented')
      end


      sgtitle(['optimization (',opt_method,') | iteration = ',num2str(iter)])  
    

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % profile modification x
      counter = counter + 1; 
      subplot(nrows, ncols, counter); hold all;
        plot(par.z/par.units.nm, x, 'r-','LineWidth',2,'DisplayName','x (no filter/window)')
        plot(par.z/par.units.nm, x_mod, 'k-','LineWidth',2,'DisplayName','x (filter + window)')
        box on
        xlabel('z (nm)')
        ylabel('x(z)')
        title('epitaxial profile modification')
        legend('Location','best')
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % full profile X = X_QW + x       
      counter = counter + 1; 
      subplot(nrows, ncols, counter); hold all;
        plot(par.z/par.units.nm, par.X_QW, 'g-','LineWidth',1,'DisplayName','X_{QW}')
        plot(par.z/par.units.nm, x_mod, 'r-','LineWidth',2,'DisplayName','x (mod)')
        plot(par.z/par.units.nm, par.X_QW + x_mod, 'b-','LineWidth',2,'DisplayName','X_{QW} + x (mod)')        
        box on
        xlabel('z (nm)')
        ylabel('X')
        title('epitaxial profile')
        legend('Location','best')
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % wave function + potential
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        U = par.U_F + par.U_QW + potential_modification(x_mod, par);

        plot(par.z/par.units.nm, U/par.units.meV , 'k-','LineWidth',2)
        plot(par.z/par.units.nm, (par.U_F + par.U_QW)/par.units.meV , 'k--','LineWidth',1)
        min_U = min(U);
        max_U = max(U);
        
        for i = 1 : par.neigs   
          if out.E(i) <= max_U
          plot(par.z/par.units.nm, out.E(i)/par.units.meV + par.N * abs(out.S(:,i)).^2, '-','Color',[1 1 1]*0.5,'LineWidth',1)
          end
        end
        i = out.idx_gnd;
        plot(par.z/par.units.nm, out.E(i)/par.units.meV + par.N * abs(out.S(:,i)).^2, '-','Color',[1 0 0],'LineWidth',2)
    
        ylim([min_U, max_U]/par.units.meV)
        box on
        xlabel('z (nm)')
        ylabel('U (meV)')
        title('wave function')
        %legend()
       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % adjoint 
      %{
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        plot(par.z/par.units.nm, par.const.e0 * out.chi , 'k-','LineWidth',2)
        box on
        xlabel('z (nm)')
        ylabel('\chi (m^{-1/2} eV^{-1})')
        title('adjoint solution')
        %legend()
      %}

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % step size alpha
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;    
        plot(1:iter, storage.alpha_hist(1:iter),'ro-')
        xlabel('iteration')
        ylabel('\alpha')
        title('step size')
        box on
        set(gca,'YScale','log')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % gradient      
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;

        scale_factor = max(abs(out.DJ))/max(abs(p));

        plot(par.z/par.units.nm, -out.DJ / scale_factor,'r-','LineWidth',2,'DisplayName','-DJ (scaled)')
        plot(par.z/par.units.nm,  p,'k-','LineWidth',2,'DisplayName','p')
        xlabel('z (nm)')
        ylabel('D_XJ(z) / p(z)')
        title('cost gradient')
        box 
        legend('Location','best')
       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % total cost
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        plot(1:iter, storage.J_tot_list(1:iter),'ko-')
        xlabel('iteration')
        ylabel('J_{tot}')
        title('cost J_{tot}')
        box on
        set(gca,'YScale','log')
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % cost J0
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        plot(1:iter, storage.J0_list(1:iter),'ko-')
        xlabel('iteration')
        ylabel('J_{0}')
        title('cost J_{0}')
        box on
        set(gca,'YScale','log')
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % cost J1 and J2       
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        plot(1:iter, storage.J1_list(1:iter),'go-','DisplayName','J_1')
        plot(1:iter, storage.J2_list(1:iter),'bo-','DisplayName','J_2')
        xlabel('iteration')
        ylabel('J_{1,2}')
        title('cost J_{1,2} (Ge-budget)')
        box on
        legend('Location','best')
        set(gca,'YScale','log')
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % cost J3      
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        plot(1:iter, storage.J3_list(1:iter),'ko-')
        xlabel('iteration')
        ylabel('J_{3}')
        title('cost J_{3} (admissible)')
        box on
        set(gca,'YScale','log')    
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % gradient norm
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        plot(1:iter, storage.DJ_sq_list(1:iter),'ko-','DisplayName','<DJ,DJ>')
        plot(1:iter, storage.p_DJ_list(1:iter),'ro-','DisplayName','-<p,DJ>')
        xlabel('iteration')
        ylabel('<DJ,DJ>')        
        title('cost gradient')
        box on
        set(gca,'YScale','log') 
        legend('Location','best')

    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % admissible domain penality function (+ derivative)
      %{
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        [f_p, df_p] = f_penalty(par.X_QW + x_mod, par);
        plot(par.z/par.units.nm, f_p,'k-','DisplayName','f')
        plot(par.z/par.units.nm, df_p,'r-','DisplayName','df')
        xlabel('z (nm)')
        ylabel('f_{penalty}')
        title('f_{penalty}')
        legend('Location','best')
        box on
        %set(gca,'YScale','log')
      %}
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % mean Ge concentration
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;    
        plot(1:iter, storage.X_mean_1(1:iter),'g-','LineWidth',2,'DisplayName','X_{mean}^{(1)}')
        plot(1:iter, storage.X_mean_2(1:iter),'b-','LineWidth',2,'DisplayName','X_{mean}^{(2)}')
        yline(par.opt.X_budget,'k--','LineWidth',2,'DisplayName','X_{budget}')
        xlabel('iteration')
        ylabel('X_{mean}')
        title('mean Ge concentration')
        box on
        legend('Location','best')
    

    
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % rician 
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        nu    = out.nu;
        sigma = out.sigma;

        Delta_det = out.Delta_det;
        Gamma     = out.Gamma;
        C         = out.C;

        mean_EVS = out.M;
        std_EVS  = sqrt(out.V);

        EVS_range = linspace(0,mean_EVS + 6*std_EVS, 101);
        pdf       = rician_pdf(EVS_range/par.units.ueV, nu/par.units.ueV, sigma/par.units.ueV);
        plot(EVS_range/par.units.ueV, pdf,'b-','LineWidth',2,'DisplayName','Rice')
        xline(mean_EVS/par.units.ueV,'r-','LineWidth',2,'DisplayName',['mean = ',num2str(out.M/par.units.ueV,'%.2f'),' ueV'])

        xline((mean_EVS + std_EVS)/par.units.ueV,'m-','LineWidth',2,'DisplayName',['mean \pm std'])
        xline((mean_EVS - std_EVS)/par.units.ueV,'m-','LineWidth',2,'DisplayName',['std = ',num2str(sqrt(out.V)/par.units.ueV,'%.2f'),' ueV'])

        xlabel('E_{VS} (ueV)')
        ylabel('pdf')
        %title(['VS stats: \nu = ',num2str(nu/par.units.ueV,'%.2f'),' ueV, \sigma = ',num2str(sigma/par.units.ueV,'%.2f'),' ueV | \nu/\sigma = ',num2str(nu/sigma,'%.3f')])
        %title(['mean(E_{VS}) = ',num2str(mean_EVS/par.units.ueV,'%.2f'),' ueV, std(E_{VS}) = ',num2str(std_EVS/par.units.ueV,'%.2f'),' ueV | Q = ',num2str(2*abs(Delta_det)/mean_EVS,'%.3f')])
        title(['Rice distribution | Q = ',num2str(2*abs(Delta_det)/mean_EVS,'%.3f'),...
          ' | |\Delta_{det}|/(2*\Gamma)^{1/2} = ',num2str(abs(Delta_det)/sqrt(2*Gamma),'%.3f')])
        box on
        legend('Location','best')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % complex plane Delta
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;    
        phi_angle = linspace(0,2*pi,101);
        
        Delta_det = out.Delta_det;
        %radius = out.sigma;
        radius = sqrt(Gamma/2);


        % plot center (complex mean value)
          plot3(real(Delta_det)/par.units.ueV, imag(Delta_det)/par.units.ueV, 1, 'ko','DisplayName',['\Delta_{det} = (',num2str(Delta_det/par.units.ueV,'%.2f'),') ueV'])

        % plot circular contour
          plot3((real(Delta_det) + sqrt(Gamma/2) * cos(phi_angle) )/par.units.ueV, ...
            (imag(Delta_det) + sqrt(Gamma/2) * sin(phi_angle) )/par.units.ueV, ...
            ones(size(phi_angle)) , ...
            'b-','LineWidth',2,'DisplayName',['circular'])
          
        % plot exact contour
          % transformation matrix
            matrix_S = 1/sqrt(2*abs(C)*(abs(C)-real(C))) * [real(C)-abs(C), imag(C); imag(C), abs(C)-real(C)];

          % main axis ellipse
            z1 = sqrt(0.5*(Gamma-abs(C))) * cos(phi_angle);
            z2 = sqrt(0.5*(Gamma+abs(C))) * sin(phi_angle);

          % result: x - Re(Delta_det) = curve(1,:), y - Im(Delta_det) = curve(2,:)
            curve = matrix_S * [z1; z2];

            plot3( (real(Delta_det) + curve(1,:))/par.units.ueV , ...
                  (imag(Delta_det) + curve(2,:))/par.units.ueV , ...
                  ones(size(phi_angle)) , ...
                  'k.','LineWidth',2,'DisplayName','exact')  

        xline(0,'k--','DisplayName','origin')  
        yline(0,'k--','DisplayName','origin')  

        legend('Location','best');

        plot_xmax = (abs(Delta_det) + 2*radius)/par.units.ueV;

        xlim([-1,1]*plot_xmax)
        ylim([-1,1]*plot_xmax)


        

        xlabel('Re(\Delta) (ueV)')
        ylabel('Im(\Delta) (ueV)')
        %title('complex plane \Delta')
        title([ ...
          ...%'\Delta_{det} = (',num2str(Delta_det/par.units.ueV,'%.2f'),') ueV |, 
          '\Gamma = ',num2str(Gamma/par.units.ueV^2,'%.2f'),' ueV^2 | ' ...
          'C = (',num2str(C/par.units.ueV^2,'%.2f'),') ueV^2' ...
          ])

        box on


  % plot pdf
    N_pts = 101;
    x = linspace(-1,+1,N_pts)*plot_xmax;
    y = linspace(-1,+1,N_pts)*plot_xmax;

    [X,Y] = meshgrid(x,y);


    pdf_mode = 2;
    switch(pdf_mode)
      case 0 % off
        p = zeros(N_pts,N_pts);
      case 1 % circular
        p = 1/(pi*Gamma/par.units.ueV^2) * exp( - 1/(Gamma/par.units.ueV^2 ) * ...
            (   (X - real(Delta_det)/par.units.ueV).^2 ...
              + (Y - imag(Delta_det)/par.units.ueV).^2 ) );
        
      case 2 % complex normal
        p = 1/(pi*sqrt(Gamma^2 - abs(C)^2)/par.units.ueV^2) * exp( - par.units.ueV^4/(Gamma^2 - abs(C)^2) * ...
            (   (Gamma-real(C))/par.units.ueV^2 *(X - real(Delta_det)/par.units.ueV).^2 ...
            + (Gamma+real(C))/par.units.ueV^2 * (Y - imag(Delta_det)/par.units.ueV).^2 ...
            - 2*imag(C)/par.units.ueV^2 * (X - real(Delta_det)/par.units.ueV) .* (Y - imag(Delta_det)/par.units.ueV) ));

      otherwise
        error('invalid option')
    end
        
    if pdf_mode > 0
    surf(X,Y,p,'DisplayName','pdf')
    shading interp

    cmap = white2red(256);
    colormap(cmap)
    end

        axis square

        
        
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Rice parameters and valley splitting moments
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;    
        plot(1:iter, storage.nu_list(1:iter)/par.units.meV,'c-','LineWidth',2,'DisplayName','\nu')
        plot(1:iter, storage.sigma_list(1:iter)/par.units.meV,'m-','LineWidth',2,'DisplayName','\sigma')
        plot(1:iter, storage.mean_EVS_list(1:iter)/par.units.meV,'r-','LineWidth',2,'DisplayName','mean')
        plot(1:iter, storage.std_EVS_list(1:iter)/par.units.meV,'b-','LineWidth',2,'DisplayName','std')
        xlabel('iteration')
        ylabel('energy (meV)')
        title('Rice parameters and VS moments')
        box on
        set(gca,'YScale','log')   
        legend('Location','best')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Fourier transforms
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        k_scale = 2*pi/par.a0;
        
        spec = abs(fft(U));
        plt_idx = (par.k>0);
        plot(par.k(plt_idx)/k_scale, spec(plt_idx)/spec(1),'k-','LineWidth',2,'DisplayName','FT(U)')
        
        spec = abs(fft(potential_modification(x_mod,par)));
        plot(par.k(plt_idx)/k_scale, spec(plt_idx)/spec(1),'b-','LineWidth',2,'DisplayName','FT(u)')
        
        spec = abs(fft(U.*out.psi0.^2));
        plot(par.k(plt_idx)/k_scale, spec(plt_idx)/spec(1),'r-','LineWidth',2,'DisplayName','FT(U \psi^2)')
        
        xline(par.opt.k_c/k_scale,'g--','LineWidth',2,'DisplayName','k_c')
        xline(2*par.k0/k_scale,'c--','LineWidth',2,'DisplayName','2 k_0')
        xline(2*(k_scale - par.k0)/k_scale,'m--','LineWidth',2,'DisplayName','2 k_1')
        
        plot(par.k(plt_idx)/k_scale, par.opt.filter(plt_idx),'g-','LineWidth',2,'DisplayName','filter ')
        ylim([1E-6,1E1])
        xlabel('k (2\pi/a_0)')
        ylabel('FT')
        title('Fourier transform')
        legend('Location','best')
        box on
        set(gca,'XScale','log')
        set(gca,'YScale','log')

      drawnow
    end

end
