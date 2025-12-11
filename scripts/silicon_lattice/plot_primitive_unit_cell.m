function [R_list] = plot_primitive_unit_cell(eps,xi)
% plot lattice
% export list of all atomic positions

  % basis vectors (unstrained)
    a{1} = 0.5*[0 1 1]';
    a{2} = 0.5*[1 0 1]';
    a{3} = 0.5*[1 1 0]';
    
  % shift between fcc lattices  
    tau = 0.25*[1 1 1]';

  % apply strain
    Id = eye(3);
    for i = 1 : 3
      a{i} = (Id + eps)*a{i};
    end   
    tau0 = (Id + eps)*tau;

   %%%%%%%%%%%%%%%%%%%%
   % list of cube vertices (1st lattice)
     n = [0 0 0;
          1 0 0;
          0 1 0;
          0 0 1];
   



    %%%%%%%%%%%%%%%%% 
    % list of edges (1st lattice, between vertices)
      e = [1 2;
           1 3;
           1 4;
           2 3;
           2 4;
           3 4;];
      
    %%%%%%%%%%%%%%%%%%%% 
    % 2nd sublattice (points in n which are shifted by tau)
      m = [1];

    % tetrahedral bonds between 1st and 2nd lattice   
      k = [1 2 3 4];





     sphere_scale = 0.1;
     [Xs,Ys,Zs] = sphere(50);

     % plot atoms of 1st lattice
       R_list = [];
     for in = 1 : size(n,1)

       R_vec = func_R1(n(in,:),a);

       R_list = [R_list, R_vec];
       
       obj = surf(R_vec(1)+sphere_scale*Xs, R_vec(2)+sphere_scale*Ys, R_vec(3)+sphere_scale*Zs);
       set(obj,'FaceColor',[1 0 0],'EdgeColor','none')

     end


     % edges of 1st fcc lattice
     for ie = 1 : size(e,1)
        Ra = func_R1(n(e(ie,1),:),a);
        Rb = func_R1(n(e(ie,2),:),a);
        [Xc Yc Zc] = cylinder2P(0.02,20,Ra', Rb');
        obj = surf(Xc,Yc,Zc);
        set(obj,'FaceColor',[1 1 1]*0.5,'EdgeColor','none')
     end
     %}

          

     % plot atoms of 2nd lattice
     for im = 1 : size(m,1)

       R_vec = func_R2(m(im),k(im,:),n,a,tau0,xi);

       R_list = [R_list, R_vec];

       obj = surf(R_vec(1)+sphere_scale*Xs, R_vec(2)+sphere_scale*Ys, R_vec(3)+sphere_scale*Zs);
       set(obj,'FaceColor',[0 0 1],'EdgeColor','none')

     end

     %
     for im = 1 : size(m,1)
        Ra = func_R2(m(im),k(im,:),n,a,tau0,xi);
        for l = 1 : 4
        Rb = func_R1(n(k(im,l),:),a);


        [Xc Yc Zc] = cylinder2P(0.02,50,Ra', Rb');
        obj = surf(Xc,Yc,Zc);
        set(obj,'FaceColor',[1 1 1]*0.5,'EdgeColor','none')
        end
     end
     %}



     axis equal
     %view([-6 14])
     view([-23 24])
     camlight
     lighting gouraud
     material shiny
     axis off
     xlim([-0.5 1.5])
     ylim([-0.5 1.5])
     zlim([-0.5 1.5])
     
end