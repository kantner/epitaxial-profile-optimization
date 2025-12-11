function [] = plot_complex_Delta(Delta_det, Gamma, C, pdf_mode, par)
% pdf_mode:
% 0 = off 
% 1 = circular
% 2 = complex normal



  % rescale to ueV
    Delta_det = Delta_det/par.units.ueV;
    Gamma     = Gamma/par.units.ueV^2;
    C         = C/par.units.ueV^2;

  % open figure
    figure(66666);clf; hold all;

    box on
    xlabel('Re(\Delta) (ueV)')
    ylabel('Im(\Delta) (ueV)')

    maxDelta = (abs(Delta_det) + 5*sqrt(Gamma));
    xlim([-1,+1]*maxDelta)
    ylim([-1,+1]*maxDelta)

    xline(0,'k--','DisplayName','origin')
    yline(0,'k--','DisplayName','origin')
    title(['intervalley coupling parameter \Delta | \Delta_{det} = (',num2str(Delta_det,'%.2f'), ') ueV, \Gamma = ',num2str(Gamma,'%.2f'),' ueV^2,  C = (',num2str(C,'%.2f'),') ueV^2'])

    axis square

    legend();

  % parameter for parametric plot  
    phi = linspace(0,2*pi,101);

    z_level = 1;

  % plot circular approximation  
    plot3( real(Delta_det) + sqrt(0.5*Gamma)*cos(phi) , ...
          imag(Delta_det) + sqrt(0.5*Gamma)*sin(phi) , ...
          z_level * ones(size(phi)), ...
          'r-','LineWidth',2,'DisplayName','circular')

  % plot true complex normal with correlation
    % We plot the contour line where the argument of the squared exponential equals 1/2.
    % This is done via diagonalization of the correlation matrix and the transformation of the resulting main axis ellipse into a transformed ellipse

    % transformation matrix
      S = 1/sqrt(2*abs(C)*(abs(C)-real(C))) * [real(C)-abs(C), imag(C); imag(C), abs(C)-real(C)];

    % main axis ellipse
      z1 = sqrt(0.5*(Gamma-abs(C))) * cos(phi);
      z2 = sqrt(0.5*(Gamma+abs(C))) * sin(phi);

    % result: x - Re(Delta_det) = curve(1,:), y - Im(Delta_det) = curve(2,:)
      curve = S * [z1; z2];

      plot3( real(Delta_det) + curve(1,:) , ...
            imag(Delta_det) + curve(2,:) , ...
            z_level * ones(size(phi)) , ...
          'b--','LineWidth',2,'DisplayName','exact')


      

  % plot pdf
    N_pts = 101;
    x = linspace(-1,+1,N_pts)*maxDelta;
    y = linspace(-1,+1,N_pts)*maxDelta;

    [X,Y] = meshgrid(x,y);


    switch(pdf_mode)
      case 0 % off
        p = zeros(N_pts,N_pts);
      case 1 % circular
        p = 1/(pi*Gamma) * exp( - 1/(Gamma) * ...
            (   (X - real(Delta_det)).^2 ...
              + (Y - imag(Delta_det)).^2 ) );
        
      case 2 % complex normal
        p = 1/(pi*sqrt(Gamma^2 - abs(C)^2)) * exp( - 1/(Gamma^2 - abs(C)^2) * ...
            (   (Gamma-real(C))*(X - real(Delta_det)).^2 ...
            + (Gamma+real(C))*(Y - imag(Delta_det)).^2 ...
            - 2*imag(C) * (X - real(Delta_det)).*(Y - imag(Delta_det)) ) );

      otherwise
        error('invalid option')
    end

    if pdf_mode > 0
    surf(X,Y,p,'DisplayName','pdf')
    shading interp

    cmap = white2red(256);
    colormap(cmap)
    end

end