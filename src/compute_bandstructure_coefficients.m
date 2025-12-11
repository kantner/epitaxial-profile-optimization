function [c, C2, C4, D4] = compute_bandstructure_coefficients(n_range, par)

    fprintf(1,'\n')
    fprintf(1,'  Compute bandstructure coefficients ...\n')
    par.eps
    %%%%%%%%%
    % plot results
      plot_results = 0;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fourier domain Bloch factors at k0
      fprintf(1,'  Fourier domain Bloch factors ... ')
      tic
      H = Hamiltonian(par.k_Delta, par);
      [c, ~] = eig(H);
      c = c(:,par.pp.idx_CB);
      toc

      par.c = c;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % tolerance
      tol = 1E-10; % tolerance for inexact evaluation of selection rule
    
    if nargout >= 2
      fprintf(1,'  coefficients C(2) .............. ')
      tic
      C2  = coeff_C2(n_range, tol, par);
      toc
    else
      fprintf(1,'  skip computation of C(2).\n')
    end

    if nargout >=3
      fprintf(1,'  coefficients C(4) .............. ')
      tic
      C4 = coeff_C4(n_range, C2);
      toc   
    else
      fprintf(1,'  skip computation of C(4).\n')
      C4 = zeros(size(n_range));
    end

    if nargout >=4
      fprintf(1,'  coefficients D(4) .............. ')
      tic
      D4 = coeff_D4(n_range, C2);
      toc   
    else
      fprintf(1,'  skip computation of D(4).\n')
      D4 = zeros(size(n_range));
    end    

    if nargout >= 2
    % print bandstructure coefficients  
      fprintf(1,'\n')
      fprintf(1,'bandstructure coefficients\n')
      switch(nargout)
        case 2
          fprintf(1,'n\tC2\n')
        case 3
          fprintf(1,'n\tC2\t\tC4\n')
        case 4
          fprintf(1,'n\tC2\t\tC4\t\tD4\n')
      end
      
      for i = 1 : length(n_range)
        switch(nargout)
          case 2
            fprintf(1,'%+d\t%+.3e\n',n_range(i), C2(i));
          case 3
            fprintf(1,'%+d\t%+.3e\t%+.3e\n',n_range(i), C2(i), C4(i));
          case 4
            fprintf(1,'%+d\t%+.3e\t%+.3e\t%+.3e\n',n_range(i), C2(i), C4(i),D4(i));
        end        
      end

    % plot  bandstructure coefficients
    if plot_results == 1      
    figure(200); clf; hold all;    
      plot(n_range, abs(C2), 'ro-','MarkerSize',10,'MarkerFaceColor','r','DisplayName','C^{(2)}')
      if nargout>=3
      plot(n_range, abs(C4), 'go-','MarkerSize',10,'MarkerFaceColor','g','DisplayName','C^{(4)}')
      end
      xlabel('n')
      ylabel('coefficient')
      title('band structure coefficients')
      box on
      legend()
      set(gca,'YScale','log')   
    end
    end

end