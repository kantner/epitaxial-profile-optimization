function [x_opt, x_mod, out] = makeFig_optimize(optimization_objective, conventional, F, X_budget, kc_list, maxiter_BB, maxiter_LBFGS, restarts, refinements, plot_skip, continuation, par)
% carry out optimization and export data + overview plot


  %%%%%%%%%%%%%%%%
  % optimization
    addpath('optimization')
    addpath('optimization/functions')
    addpath('optimization/user')

    

  % cost functional J0: weights and parameters
    % optimization_objective ... % 1 = min std/mean | 2 = max mean | 3 = min std | 4 = max std | 5 = min mean
    [par] = set_J0_parameters(optimization_objective, par);



  % electric field (in V/m)
    par.F    = F;
    par.U_F  = -par.const.e0 * par.F * par.z;    
        

  % continuation_mode
    continuation_mode = continuation; % 0 = start from zero | 1 = from previous




    

for i_kc = 1 : length(kc_list)
  if conventional == 1 % use standard cosine profile

    % set some dummy values because we need to
      par.opt.window       = par.QW_indicator;
      par.opt.k_c          = 10*2*pi/par.a0;%
      par.opt.filter_type  = 0; % 0 = no filter | 1 = rect | 2 = exp | 3 = gauss | 4 = bessel | 5 = fermi | 6 = erf
      par.opt.apply_filter = 0;
      par.opt.apply_window = 0;
      par.opt.filter       = generate_filter(par.k, par.opt.k_c, par.opt.filter_type);

      par.opt.w = [0 0 0];
      
      %h = par.h_QW;
      q = kc_list(i_kc);
      %X_WW = X_budget / (0.5*(1+sin(q*h)/(q*h)));
      %X_WW = X_budget 
      %q = 2*par.k1;
      phi = 0;
      [X_WW] = compute_wiggle_well_amplitude(X_budget, q,phi,par);

      x = x_wiggle_well(X_WW, q, 0, par);

      x_mod = x;

    % set dummy output
      x_opt = x;

      

    
      




    else % optimization


  %  par.opt.a = [0 0 1 0];
  %  par.opt.E = [1 1]*par.units.meV;


  % window
    par.opt.apply_window = 1; % 0 = off | 1 = on    
    par.opt.window       = par.QW_indicator;

  % filter
    par.opt.apply_filter = 1; % 0 = off | 1 = on
    par.opt.filter_type  = 1; % 0 = no filter | 1 = rect | 2 = exp | 3 = gauss | 4 = bessel | 5 = fermi | 6 = erf
    par.opt.k_c          = kc_list(i_kc); %
    par.opt.filter       = generate_filter(par.k, par.opt.k_c, par.opt.filter_type);

  % Ge budget
    par.opt.X_budget = X_budget;

  % penalty function
    par.opt.penalty.eps   = 0 * 0.001;
    par.opt.penalty.x_max = 0.3 + 1*par.opt.penalty.eps;
    par.opt.penalty.x_min = 0.0 - 1*par.opt.penalty.eps;

  % plot penalty function
    %{
    x_range = linspace(par.opt.penalty.x_min - abs(par.opt.penalty.x_max-par.opt.penalty.x_min)/par.opt.penalty.eps,par.opt.penalty.x_max+abs(par.opt.penalty.x_max-par.opt.penalty.x_min)/par.opt.penalty.eps,100001);
    figure(13243); clf; hold all;
    [f, df] = f_penalty(x_range, par);
    plot(x_range,f,'bo-')
    plot(x_range,df,'ro-')
    xline(x_range(end),'c-')
    xline(x_range(1),'c-')
    %}

  % dummy value for invalid result (negative <|Delta_rand|^2> for invalid X)
    par.sigma_failure = 1E100;

  % print info  
    par.opt.print_info = 2;

  % numerics
    par.neigs = 10;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % optimization method
    par.opt.method     = 2;  % 1 = gradient descend | 2 = Barzilai-Borwein | 3 = L-BFGS
    par.opt.maxiter    = 5; 
    par.opt.f_tol      = -Inf;
    par.opt.dfnorm_tol = -Inf;
    par.opt.dx_upd_tol = -Inf; % not good    
    par.opt.df_upd_tol = -Inf; % not good
    par.opt.f_upd_tol  = 1E-10; % detect stagnation
    par.opt.n_stagnate = 10; % detect stagnation
    par.opt.print_info = 2; % 0 = off | 1 = basic info | 2 = detailed info
    par.opt.plot_skip  = 1; % update plot after plot_skip steps 
    par.opt.adjoint_method = 4; % method for adjoint solution (1..3: Green's function | 4 = backslash + orthogonality projection | 5 = pseudoinverse | 6 = augmentation)
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % method specific parameters

  % Armijo linesearch parameters
    par.opt.Armijo.alpha0  = 1E-13;
    par.opt.Armijo.maxiter = 20;
    par.opt.Armijo.c1      = 0.5;
    par.opt.Armijo.shrink  = 0.5;
    par.opt.Armijo.increase = 1.5;
    par.opt.Armijo.print_info = 2; % 0 = off | 1 = on
    par.opt.Armijo.plot_linesearch  = 0; % 0 = off | 1 = on
    par.opt.Armijo.plot_alpha_range = logspace(-40,-16,21);
    par.opt.Armijo.plot_alpha_scale = 'log';  % 'linear' or 'log'
    
  % Wolfe linesearch parameters
    par.opt.Wolfe.alpha_init = 1E-11; % only first step
    par.opt.Wolfe.alpha0     = 1;     % standard for scaled updates
    par.opt.Wolfe.maxiter = 20;
    par.opt.Wolfe.c1      = 0.5;
    par.opt.Wolfe.c2      = 0.9;
    par.opt.Wolfe.lambda  = 2;
    par.opt.Wolfe.print_info = 2; % 0 = off | 1 = on    
  
  % gradient descent parameters
    par.opt.gradientDescent.linesearch_mode = 1; % 1 = Armijo | 2 = Wolfe
    
  % Barzilai Borwein parameters    
    par.opt.BB.alpha0     = 1E-13;
    par.opt.BB.dx_max     = 1E-3;
    par.opt.BB.alpha_mode = 6; % mode for step size computation: 1 = long | 2 = short | 3 = geom. mean | 4 = alternating | 5 = random | 6 = positive alternating
    
  % L-BFGS parameters
    %par.opt.LBFGS.memory   = 1000;
    %par.opt.LBFGS.memory   = 20;
    par.opt.LBFGS.max_skip = 25;
    par.opt.LBFGS.min_curv = 1E-19;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize epitaxial profile  
    if continuation_mode == 0
      x = zeros(par.N, 1);


      %x = exp( -0.5*((par.z+4*par.units.nm)/(1*par.units.nm)).^2 );

      
    % initialize constant with target Ge budget
      x = par.QW_indicator;
      x = sum(par.QW_indicator) *  X_budget/sum(x) * x;


    %x = 0 * 0.01 * 0.5*(1+cos(2*(2*pi/par.a0 - par.k0) *  par.z )) .* par.QW_indicator;
    else
      if i_kc == 1
        x = zeros(par.N, 1);
      else
        x = x_opt;
      end
    end

  % measure performance
  %  profile clear
  %  profile on
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % optimization   
  
  % penalty weights  
    par.opt.w = [0 1E4 1E5];   % w1 (Ge-Budget), w2 (Ge-Budget), w3 (max/min Ge limits)

  % stage 1: BB
    if maxiter_BB > 0
    par.opt.method    = 2; % 1 = gradient descent | 2 = Barzilai-Borwein | 3 =  L-BFGS
    par.opt.maxiter   = maxiter_BB; 
    par.opt.plot_skip = plot_skip;
    [x_opt, stats, storage] = optimization(x, par, par.opt);

    x = x_opt;
    end

  % stage 2: LBFGS
  %{
    par.opt.LBFGS.memory   = 10;
    par.opt.method    = 3; % 1 = gradient descent | 2 = Barzilai-Borwein | 3 =  L-BFGS
    par.opt.maxiter   = maxiter_LBFGS; 
    par.opt.plot_skip = plot_skip;
    [x_opt, stats, storage] = optimization(x, par, par.opt);

    x = x_opt;
    x = max(x + randn(size(x)),0);
%}
  % stage 2: LBFGS
    for j = 1 : restarts
    par.opt.LBFGS.memory   = 20;
    par.opt.method    = 3; % 1 = gradient descent | 2 = Barzilai-Borwein | 3 =  L-BFGS
    par.opt.maxiter   = maxiter_LBFGS; 
    par.opt.plot_skip = plot_skip;
    [x_opt, stats, storage] = optimization(x, par, par.opt);

    x = x_opt;
    %x = max(x + 0.01*randn(size(x)),0);    
    end
    
  % refinement stage: LBFGS + increased penalties
    for j = 1 : refinements
    par.opt.LBFGS.memory   = maxiter_LBFGS;
    par.opt.w(2)      = 10 * par.opt.w(2);   % w1 (Ge-Budget), w2 (Ge-Budget), w3 (max/min Ge limits)
    par.opt.method    = 3; % 1 = gradient descent | 2 = Barzilai-Borwein | 3 =  L-BFGS
    par.opt.maxiter   = maxiter_LBFGS; 
    par.opt.plot_skip = plot_skip;
    [x_opt, stats, storage] = optimization(x, par, par.opt);

    x = x_opt;

    end

    
      end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % make plot
    fig_obj = figure(9999); clf; hold all
      nrows = 2;
      ncols = 3;

      sgtitle(['results (opt-obj-id = ',num2str(optimization_objective),' | Ge budget = ',num2str(X_budget*100),' % | k_c = ',num2str(par.opt.k_c*par.a0/(2*pi)),' \times 2\pi/a_0)'])


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

      % solve again at optimum to provide out-struct
        [~, out]  = function_f(x, par);
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % full profile X = X_QW + x       
      counter = 1;
      subplot(nrows, ncols, counter); hold all;
        plot(par.z/par.units.nm, par.X_QW, '--','Color',[1 1 1]*0.5,'LineWidth',1,'DisplayName','X_{QW}')
        %plot(par.z/par.units.nm, x_mod, 'r-','LineWidth',2,'DisplayName','x (mod)')
        plot(par.z/par.units.nm, par.X_QW + x_mod, 'b-','LineWidth',2,'DisplayName','X_{QW} + x (mod)')        
        box on
        xlabel('z (nm)')
        ylabel('X')
        title('epitaxial profile')
        legend('Location','best')
        %xlim([-15.5,5.5])
        %ylim([-0.025,0.375])
        xlim([-11,11])
        ylim([-0.025,0.375])


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
        %xlim([-11,11])        
        %ylim([-35, 215])
        xlim([-11,11])
        ylim([-35, 215])

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Rician 
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        nu    = out.nu;
        sigma = out.sigma; % sigma = sqrt(2*Gamma)

        Delta_det = out.Delta_det;
        Gamma     = out.Gamma;
        C         = out.C;

        mean_EVS = out.M;
        std_EVS  = sqrt(out.V);


        zeta = abs(Delta_det)/sqrt(2*Gamma);

        Q = nu/mean_EVS;

       
        if (sigma - sqrt(2*Gamma)) > 1E-12
          error('sigma and Gamma do not match!')
        end

        %EVS_range = linspace(0,mean_EVS + 6*std_EVS, 101);
        EVS_range = linspace(0, 3.1, 1001)*par.units.meV;
        pdf       = rician_pdf(EVS_range/par.units.ueV, nu/par.units.ueV, sigma/par.units.ueV);
        plot(EVS_range/par.units.ueV, pdf,'b-','LineWidth',2,'DisplayName','Rice')
        xline(mean_EVS/par.units.ueV,'r-','LineWidth',2,'DisplayName',['mean = ',num2str(out.M/par.units.ueV,'%.2f'),' ueV'])

        xline((mean_EVS + std_EVS)/par.units.ueV,'m-','LineWidth',2,'DisplayName',['mean \pm std'])
        xline((mean_EVS - std_EVS)/par.units.ueV,'m-','LineWidth',2,'DisplayName',['std = ',num2str(sqrt(out.V)/par.units.ueV,'%.2f'),' ueV)'])

        xlabel('E_{VS} (ueV)')
        ylabel('pdf')
        title(['\nu = ',num2str(nu/par.units.ueV,'%.2f'),' ueV, \sigma = ',num2str(sigma/par.units.ueV,'%.2f'),' ueV | \zeta = ',num2str(zeta,'%.3f'),' | Q = ',num2str(Q,'%.3g')])
        box on
        legend('Location','best')        
        
        %xlim([-0.1, 3.1]*1E3)
        %ylim([-0.05, 2.15]*1E-3)
        xlim([-0.1, 1.1]*1E3)
        ylim([-0.05, 7.05]*1E-3)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % complex plane Delta
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;    
        phi_angle = linspace(0,2*pi,101);
        
        
        plot_xmax = 450;

%{
        plot(real(Delta_det)/par.units.ueV, imag(Delta_det)/par.units.ueV, 'ro')
        plot((real(Delta_det) + radius * cos(phi_angle) )/par.units.ueV, (imag(Delta_det) + radius * sin(phi_angle) )/par.units.ueV, 'r--','LineWidth',2)

        xline(0,'k--')
        yline(0,'k--')

        plot_xmax = 1100;

        xlim([-1,1]*plot_xmax)
        ylim([-1,1]*plot_xmax)

        xlabel('Re(\Delta) (ueV)')
        ylabel('Im(\Delta) (ueV)')
        title('complex plane \Delta')
        box on    
        axis square
%}


        % plot center (complex mean value)
          plot3(real(Delta_det)/par.units.ueV, imag(Delta_det)/par.units.ueV, 1, 'k.','DisplayName',['\Delta_{det} = (',num2str(Delta_det/par.units.ueV,'%.2f'),') ueV'])

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
                  'k-.','LineWidth',2,'DisplayName','exact')  

        xline(0,'k--','DisplayName','origin')  
        yline(0,'k--','DisplayName','origin')  

        legend('Location','best');
       

        xlim([-1,1]*plot_xmax)
        ylim([-1,1]*plot_xmax)


        

        xlabel('Re(\Delta) (ueV)')
        ylabel('Im(\Delta) (ueV)')
        %title('complex plane \Delta')
        title([ ...
          ...%'\Delta_{det} = (',num2str(Delta_det/par.units.ueV,'%.2f'),') ueV |, 
          '\Gamma^{(1/2)} = ',num2str(sqrt(Gamma/par.units.ueV^2),'%.2f'),' ueV | ' ...
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
      % Fourier transforms
      counter = counter + 1;  
      subplot(nrows, ncols, counter); hold all;
        k_scale = 2*pi/par.a0;
        
        %spec = abs(fft(U));
        plt_idx = (par.k>0);
        %plot(par.k(plt_idx)/k_scale, spec(plt_idx)/spec(1),'k-','LineWidth',2,'DisplayName','FT(U)')
        
        spec = abs(fft(potential_modification(x_mod,par)));
        plot(par.k(plt_idx)/k_scale, spec(plt_idx)/spec(1),'b-','LineWidth',2,'DisplayName','FT(u)')
        
        spec = abs(fft(U.*out.psi0.^2));
        plot(par.k(plt_idx)/k_scale, spec(plt_idx)/spec(1),'r-','LineWidth',2,'DisplayName','FT(U \psi^2)')
        
        xline(par.opt.k_c/k_scale,'g--','LineWidth',2,'DisplayName','k_c')
        xline(2*par.k0/k_scale,'c--','LineWidth',2,'DisplayName','2 k_0')
        xline(2*par.k1/k_scale,'m--','LineWidth',2,'DisplayName','2 k_1')
        
        plot(par.k(plt_idx)/k_scale, par.opt.filter(plt_idx),'g-','LineWidth',2,'DisplayName','filter ')        
        xlabel('k (2\pi/a_0)')
        ylabel('FT')
        title('Fourier transform')
        legend('Location','best')
        box on
        set(gca,'XScale','log')
        set(gca,'YScale','log')

        ylim([5E-5, 5E0])
        xlim([0.8E-2, 2.5])




      drawnow

    %%%%%%%%%%%%%
    % save
      height = 800;
      width  = 1600;
      set(fig_obj,'Position',[100 100 width height]);

      %file_name = ['profile_obj_',num2str(optimization_objective),'_kc_',num2str(i_kc),'.pdf'];
      %file_name = ['profile_obj_',num2str(optimization_objective),'_kc_',num2str(kc_list(i_kc)*par.a0/(2*pi),'%.4f'),'.pdf'];
      if conventional == 0
        file_name = ['profile_obj_',num2str(optimization_objective),'_X_',num2str(X_budget,'%.2f'),'_kc_',num2str(kc_list(i_kc)*par.a0/(2*pi),'%.4f'),'.pdf'];
      elseif conventional == 1
        file_name = ['profile_obj_',num2str(optimization_objective),'_X_',num2str(X_budget,'%.2f'),'_kc_',num2str(kc_list(i_kc)*par.a0/(2*pi),'%.4f'),'_conventional.pdf'];
      end
       
      exportgraphics(fig_obj, file_name, 'ContentType','vector')




end
