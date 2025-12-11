function [R_list] = plot_diamond_lattice(eps,xi)
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
          -1 +1 +1;
          0 0 2;
          +1 -1 +1;
          +1 +1 -1;
          0 2 0;
          +1 +1 +1;
          2 0 0];
   % add face centers
     n = [n;
          0 0 1;
          0 1 0;
          0 1 1;
          1 0 1;
          1 0 0
          1 1 0;
          ];


    %%%%%%%%%%%%%%%%% 
    % list of edges (1st lattice, between vertices)
      e = [1 2;
           2 3;
           3 4;
           4 1;
           1 5;
           2 6;
           3 7;
           4 8;
           5 6;
           6 7;
           7 8;
           8, 5];
      % add edges between face centers
      %{
        e = [e;
             9 1;
             9 2;
             9 3;
             9 4;
             10 5;
             10 6;
             10 2;
             10 1;
             11 2;
             11 3;
             11 7;
             11 6;
             12 3;
             12 4;
             12 8;
             12 7;
             13 1;
             13 4;
             13 8;
             13 5;
             14 5;
             14 6;
             14 7;
             14 8];
      %}

     
    %%%%%%%%%%%%%%%%%%%% 
    % 2nd sublattice (points in n which are shifted by tau)
      m = [1; 9; 10; 13];      

    % tetrahedral bonds between 1st and 2nd lattice   
      k = [1 9 13 10;
           9 3 12 11;
           10 11 6 14;
           13 12 8 14];

      % edge list: face centers 



     sphere_scale = 0.1;
     [Xs,Ys,Zs] = sphere(20);

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


        [Xc Yc Zc] = cylinder2P(0.02,20,Ra', Rb');
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