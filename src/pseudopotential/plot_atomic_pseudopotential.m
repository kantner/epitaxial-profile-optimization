function [] = plot_atomic_pseudopotential(par)
    
    q_range = linspace(0,7,1001) * 2*pi/par.a0;
    V_loc = zeros(size(q_range));
    for iq = 1 : length(q_range)
      V_loc(iq) = local_pseudopotential(q_range(iq), par);
    end
   
    figure(100);clf;hold all;
      plot(q_range*par.a0/(2*pi), V_loc/par.units.Ry, 'ko-')
      xline(sqrt(3),'r-')
      xline(sqrt(8),'g-')
      xline(sqrt(11),'b-')
      xline(3*par.pp.kF*par.a0/(2*pi),'m-')
      yline(par.pp.V_sqrt3/par.units.Ry,'r-')
      yline(par.pp.V_sqrt8/par.units.Ry,'g-')
      yline(par.pp.V_sqrt11/par.units.Ry,'b-')
      yline(par.pp.V_0/par.units.Ry,'c-')
      yline(0,'m-')

      xlabel('q (2\pi/a_0)')
      ylabel('E (Ry)')
      box on
      drawnow

end
