clear all
clc

%% lattice vectors
   a{1} = 0.5*[0 1 1]';
   a{2} = 0.5*[1 0 1]';
   a{3} = 0.5*[1 1 0]';
   tau = 0.25*[1 1 1]';

 % strain
   eps = zeros(3,3);
   %eps(1,1) = 0.1;
   %eps(2,2) = eps(1,1);
   %eps(3,3) = -0.1;

   Id = eye(3);
 % modify lattice vectors
   for i = 1 : 3
     a{i} = (Id + eps)*a{i};
   end

   figure(1);clf;hold all;
   
   % cube vertices
     n = [0 0 0;
          -1 +1 +1;
          0 0 2;
          +1 -1 +1;
          +1 +1 -1;
          0 2 0;
          +1 +1 +1;
          2 0 0];

   % face centers
     %
     n = [n;
          0 0 1;
          0 1 0;
          0 1 1;
          1 0 1;
          1 0 0
          1 1 0;
          ];
     %}

    % edge list: vertices
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

      % edge list: face centers
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

      R = @(n) n(1)*a{1} + n(2)*a{2} + n(3)*a{3};

     sphere_scale = 0.1;
     [X,Y,Z] = sphere(20);
     for in = 1 : size(n,1)

       %R = n(in,1)*a{1} + n(in,2)*a{2} + n(in,3)*a{3};
       R_vec = R(n(in,:));
       
       obj = surf(R_vec(1)+sphere_scale*X, R_vec(2)+sphere_scale*Y, R_vec(3)+sphere_scale*Z);
       set(obj,'FaceColor',[1 1 1]*0.5,'EdgeColor','none')

     end

     for ie = 1 : size(e,1)
        R1 = R(n(e(ie,1),:));
        R2 = R(n(e(ie,2),:));
        %plot3([R1(1) R2(1)], [R1(2) R2(2)], [R1(3) R2(3)], 'r-','LineWidth',12)
        [X Y Z] = cylinder2P(0.02,20,R1', R2');
        obj = surf(X,Y,Z);
        set(obj,'FaceColor',[1 0 0]*0.99,'EdgeColor','none')
     end

%[X,Y,Z] = cylinder(0.025);
%surf(X,Y,Z)

     axis equal
     view(3)
     camlight
     lighting gouraud
     material shiny