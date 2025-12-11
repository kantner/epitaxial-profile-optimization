function [] = plot_bandstructure(par)

   n_steps = 101;

   n_energy = par.pp.N_G;
   if par.pp.SOI == 1
    n_energy = 2 * n_energy;
   end

   color_CB = [1 0 0];
   color_HH = [0 0 1];
   color_LH = [0 1 1];
   color_SO = [1 0 1];

%%%%%%%%%%%%%%%   
   s = linspace(0,1,n_steps);
   energy = zeros(n_energy, n_steps);




   % list of high symmetry points (relaxed)
   %{
    k_L     = [0.5 0.5 0.5]' * 2*pi/par.a0;
    k_Gamma = [0 0 0]'       * 2*pi/par.a0;
    k_X     = [0 0 1]'       * 2*pi/par.a0;
    k_K     = [3/4, 3/4, 0]' * 2*pi/par.a0;        
    k_U     = [1/4, 1/4, 1]' * 2*pi/par.a0;
    k_W     = [1/2, 0, 1]'   * 2*pi/par.a0;
   %}
 
 % path definition   
   k_start_list = [par.k_L, par.k_Gamma, par.k_X, par.k_K];
   k_final_list = [par.k_Gamma, par.k_X, par.k_U, par.k_Gamma];

   n_paths = size(k_start_list,2);

   path_length= zeros(n_paths, 1);

   for np = 1 : n_paths
     path_length(np) = sqrt(sum((k_final_list(:,np) - k_start_list(:,np)).^2)) * par.a0/(2*pi);
   end
    
  %%%%%%% 
  % top valence band energy at Gamma point (as reference energy)
    k = par.k_Gamma;
    H = Hamiltonian(k, par);
    [~,E] = eig(H);
    E = diag(E) * par.energy_scale;
    Ev_offset = E(par.pp.idx_HH);

      
   figure(101);clf;hold all;
    ylim([-12.5, 6.5])
    xlim([0, sum(path_length)])
    ylabel('energy (eV)')
    xlabel('wave vector')
    set(gca,'XTick',[0;cumsum(path_length)]')
    set(gca,'XTickLabel',{'L','\Gamma','X','K;U','\Gamma'})
    title('Si bandstructure')
    box on

   for np = 1 : n_paths
     k_start = k_start_list(:,np);
     k_final = k_final_list(:,np);     

  
     for i = 1 : n_steps
     
     si   = s(i);
     k   = (1-si) *  k_start + si * k_final;
  
     H = Hamiltonian(k, par);
  
   % diagonalization
     [~,E] = eig(H);
     E = diag(E) * par.energy_scale;
  
     energy(:,i) = E;
     end
  
     for l = 1 : n_energy
  
       switch(l)
         case par.pp.idx_CB
           color = color_CB;
         case par.pp.idx_HH
           color = color_HH;
         case par.pp.idx_LH
           color = color_LH;
         case par.pp.idx_SO
           color = color_SO;         
         otherwise
           color = [0 0 0];
       end
         
      plot(sum(path_length(1:np-1)) + s*path_length(np),(energy(l,:)-Ev_offset)/par.units.eV,'-','Color',color,'LineWidth',2)
     end
  
   % add vertical line
     if np < n_paths
     xline(sum(path_length(1:np)),'k--')
     end
  
  
     
   end
   
   drawnow

end