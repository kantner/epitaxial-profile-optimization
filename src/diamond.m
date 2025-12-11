function [a,b,tau] = diamond(a0, eps, xi)
% INPUT
%   a0  ... lattice constant
%   eps ... 3x3 strain tensor
%   xi  ... Kleinman's internal strain parameter
% OUTPUT
%   a ... lattice vectors
%   b ... reciprocal lattice vectors

  % plot option
    plot_switch = 0; % 0 = off | 1 = on 
    
  % print bond length information 
    print_bond_length = 0;



 %% no strain
  % lattice vectors
    a{1} = a0 * 0.5 * [0 1 1]';
    a{2} = a0 * 0.5 * [1 0 1]';
    a{3} = a0 * 0.5 * [1 1 0]';

  % relative shift vector between fcc lattices  
    tau = a0 * 0.25 * [1 1 1]';

  % normalizing volume (scalar triple product)
  % remark: volume of primitive cell is vol_prim = vol/(2*3)
    vol = a{1}' * cross(a{2}, a{3});

  % reciprocal lattice vectors
    b{1} = 2*pi/vol * cross(a{2},a{3});
    b{2} = 2*pi/vol * cross(a{3},a{1});
    b{3} = 2*pi/vol * cross(a{1},a{2});


  %%%%%%%%%%%%%%%%%%%%%%%%%%  
    a{4} = [0 0 0]'; % 4th vertex atom at origin

  % plot primitive cell
    if plot_switch == 1    
    figure(888); clf; hold all;    
    % vertex atoms
      color = [0 0 0];
      for i = 1 : 4
        plot3(a{i}(1)/a0, a{i}(2)/a0, a{i}(3)/a0, 'o','Color',color,'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',color)
        for j = i+1 : 4
          plot3([a{i}(1) a{j}(1)]/a0, [a{i}(2) a{j}(2)]/a0, [a{i}(3) a{j}(3)]/a0, '-','Color',color,'LineWidth',2)
        end
      end
    % central atoms
      color = [1 0 0];
      plot3(tau(1)/a0, tau(2)/a0, tau(3)/a0, 'o','Color',color,'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',color)
      for i = 1 : 4
        plot3([a{i}(1) tau(1)]/a0, [a{i}(2) tau(2)]/a0, [a{i}(3) tau(3)]/a0, '-','Color',color,'LineWidth',2)
      end       
      view(3);
      xlabel('x (a_0)')
      ylabel('y (a_0)')
      zlabel('z (a_0)')
      axis equal
      box on

    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%    
  % echo bond length
    if print_bond_length == 1
      fprintf(1,'\nbond length (relaxed)\n')
      for i = 1 : 4
        for j = i+1 : 4
        fprintf(1,'%d-%d\t%.4f\n',i,j,sqrt((a{i}-a{j})'*(a{i}-a{j}))/a0);
        end
      end
      for i = 1 : 4
      fprintf(1,'%d-%d\t%.4f\n',i,5,sqrt((a{i}-tau)'*(a{i}-tau))/a0);
      end
    end
   
 %% include strain
  % strained lattice vectors
    Id = eye(3);
    for i = 1 : 3
      a{i} = (Id + eps)*a{i};
    end
   
  % strained volume  
    vol = a{1}' * cross(a{2}, a{3});

  % strained reciprocal lattice vectors
    b{1} = 2*pi/vol * cross(a{2},a{3});
    b{2} = 2*pi/vol * cross(a{3},a{1});
    b{3} = 2*pi/vol * cross(a{1},a{2});

  %%%%%%%%%%%%%%%%%%%%%%%%%%  
  % find position of central atom
  % (A) midpoint according to macroscopic strain
    tau_A = 0.25*(a{1}+a{2}+a{3}+a{4});

  % (B) point with equal distance to vertex atoms
    mat = [(a{2}-a{1})';
           (a{3}-a{1})';
           (a{4}-a{1})'];
    
    rhs = 0.5 * [a{2}'*a{2} - a{1}'*a{1};
                 a{3}'*a{3} - a{1}'*a{1};
                 a{4}'*a{4} - a{1}'*a{1}];

    tau_B = mat\rhs;

  % convex combination using internal strain parameter
    tau = (1-xi) * tau_A + xi * tau_B;


%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % plot primitive cell
    if plot_switch == 1
    figure(888); hold all;    
    % vertex atoms
      color = [0 0 1];
      for i = 1 : 4
        plot3(a{i}(1)/a0, a{i}(2)/a0, a{i}(3)/a0, 'o','Color',color,'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',color)
        for j = i+1 : 4
          plot3([a{i}(1) a{j}(1)]/a0, [a{i}(2) a{j}(2)]/a0, [a{i}(3) a{j}(3)]/a0, '-','Color',color,'LineWidth',2)
        end
      end
    % central atoms
      color = [1 0 1];
      plot3(tau(1)/a0, tau(2)/a0, tau(3)/a0, 'o','Color',color,'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',color)
      for i = 1 : 4
        plot3([a{i}(1) tau(1)]/a0, [a{i}(2) tau(2)]/a0, [a{i}(3) tau(3)]/a0, '-','Color',color,'LineWidth',2)
      end       
      view(3);
      xlabel('x (a_0)')
      ylabel('y (a_0)')
      zlabel('z (a_0)')
      axis equal
      box on
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%    
  % echo bond lengthc
    if print_bond_length == 1
      fprintf(1,'\nbond length (strained)\n')
      for i = 1 : 4
        for j = i+1 : 4
        fprintf(1,'%d-%d\t%.4f\n',i,j,sqrt((a{i}-a{j})'*(a{i}-a{j}))/a0);
        end
      end
      for i = 1 : 4
      fprintf(1,'%d-%d\t%.4f\n',i,5,sqrt((a{i}-tau)'*(a{i}-tau))/a0);
      end  
    end


  %%%%%%%%%%%%%%%%%%%%%%%%%
  % output
  % delete atom at origin
    a(4) = [];



end