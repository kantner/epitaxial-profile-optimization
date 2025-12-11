function [par] = reciprocal_lattice_vectors(par)


%% generate list of G-vectors   
   G_list   = [];
   G_max_sq = (2*par.const.m0)/par.const.hbar^2 * par.pp.E_cutoff;
   g_max_sq = G_max_sq * (par.a0/(2*pi))^2;
   
 % set max integer for search sweep
   n_max = ceil(sqrt(G_max_sq)*(par.a0/(2*pi)));

 % rescale
   for i = 1 : 3
     b{i} = par.pp.b{i}*par.a0/(2*pi);
   end


   for nx = -n_max : n_max  
     g1 = nx * b{1};
     for ny = -n_max : n_max
       g2 = g1 + ny * b{2};
       for nz = -n_max : n_max
         g = g2 + nz * b{3};

       % check cutoff criterion and append to list       
         if g'*g < g_max_sq
           G_list = [G_list, g * 2*pi/par.a0];
         end

       end
     end
   end

   par.pp.G_list = G_list;

   par.pp.N_G = size(par.pp.G_list,2);

end