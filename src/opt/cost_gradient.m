function [out] = cost_gradient(x, out, par)
% functional derivative of cost functional w.r.t. X
% INPUT
%   x ..... epitaxial profile modification
%   out ... output struct from cost_functional (which contains functional derivatives)

  DEBUG = 0;

  %%%%%%%%%%%%%%%%%%%%%
  % apply filter and window
  % Replace x --> window * conv(filter, x) at the beginning of the
  % functional computation. This modification is passed to all subroutines.

    if isfield(par,'opt')
    % apply filter
      if par.opt.apply_filter == 1
        [x] = apply_filter(x, par.opt.filter);
      end
  
    % apply window
      if par.opt.apply_window == 1
        [x] = apply_window(x, par.opt.window);
      end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % contribution from PDE constraint

  % solve adjoint problem
    chi = adjoint_equation(x, out, par);

  % DEBUG  
    DEBUG = 0;
    if DEBUG == 1
    par.opt.adjoint_method = 1;
    tic
    chi1 = adjoint_equation(x, out, par);
    toc

    par.opt.adjoint_method = 2;
    tic
    chi2 = adjoint_equation(x, out, par);
    toc

    par.opt.adjoint_method = 3;
    tic
    chi3 = adjoint_equation(x, out, par);
    toc

    par.opt.adjoint_method = 4;
    tic
    chi4 = adjoint_equation(x, out, par);
    toc

    par.opt.adjoint_method = 5;
    tic
    chi5 = adjoint_equation(x, out, par);
    toc
  
    par.opt.adjoint_method = 6;
    tic
    chi6 = adjoint_equation(x, out, par);
    toc    


    chi_compare = [chi1, chi2, chi3, chi4, chi5, chi6]
    difference  = [ chi2-chi1, chi3-chi1, chi4-chi1, chi5-chi1, chi6-chi1]

    pause
    end


  % gradient w.r.t. X
    DJ = par.dEc * out.psi0 .* chi;% / par.dz;

  % plot
    if DEBUG == 1
      figure(11111111);clf;hold all;
        title('gradient of J w.r.t x')
        box on
        legend()
        xlabel('z (nm)')
  
        plot(par.z/par.units.nm, DJ, 'bo-','DisplayName','D_xJ_\chi')
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % contribution from J0
    if sum(abs(par.opt.a) + abs(par.opt.b))>0
      %{
      dJ0_dx =   (out.dJ0_dM * out.dM_dnu    + out.dJ0_dV * out.dV_dnu)    * out.D_nu_dX ...
               + (out.dJ0_dM * out.dM_dsigma + out.dJ0_dV * out.dV_dsigma) * out.D_sigma_dX;
      %}
      dJ0_dx =    out.dJ0_dnu    * out.D_nu_dX ...
                + out.dJ0_dsigma * out.D_sigma_dX;

      DJ = DJ + dJ0_dx;

      % plot
      if DEBUG == 1
        plot(par.z/par.units.nm, dJ0_dx, 'ko-','DisplayName','D_xJ_0')
      end

    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % contribution from J1
    if par.opt.w(1) ~= 0
      DJ = DJ + out.dJ1_dx;

      % plot
      if DEBUG == 1      
        plot(par.z/par.units.nm, out.dJ1_dx, 'ro-','DisplayName','D_xJ_1')
      end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % contribution from J2  
    if par.opt.w(2) ~= 0
      DJ = DJ + out.dJ2_dx;

      % plot
      if DEBUG == 1      
        plot(par.z/par.units.nm, out.dJ2_dx, 'go-','DisplayName','D_xJ_2')
      end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % contribution from J3
    if par.opt.w(3) ~= 0
      DJ = DJ + out.dJ3_dx;

      % plot
      if DEBUG == 1      
        plot(par.z/par.units.nm, out.dJ3_dx, 'mx-','DisplayName','D_xJ_3')
      end
    end
   

  %%%%%%%%%%%%%%%%%%%%%
  % apply filter and window to gradient (reverse order)
  % Replace DJ --> conv(filter, window * DJ) at the end of the functional
  % gradient computation to obtain gradient for the original variable.

    if isfield(par,'opt')
    % apply window
      if par.opt.apply_window == 1
        [DJ] = apply_window(DJ, par.opt.window);
      end
      
    % apply filter
      if par.opt.apply_filter == 1
        [DJ] = apply_filter(DJ, par.opt.filter);
      end
    end


  %%%%%%%%%%%%%%
  % add to out
    out.DJ  = DJ;
    out.chi = chi;
   


end