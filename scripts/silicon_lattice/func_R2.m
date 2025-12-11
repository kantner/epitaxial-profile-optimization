function R = func_R2(m,k,n,a,tau0,xi)
% find position on 2nd fcc lattice
% m ... index of original point in 1st lattice that is shifted by tau
% k ... indices of neighbors in 1st lattice
% n ... weights to find 1st lattice vectors
% a ... lattice vectors (strained)
% tau0 ... macroscopic strained tau
% xi ..... internal strain parameter

% step 1: original point with macroscopic shift
  tau0 = func_R1(n(m,:),a) + tau0;

% obtain positions of all neighbors
  for i = 1 : 4
    X{i} = func_R1(n(k(i),:),a);
  end
  %{
  X{1}
  X{2}
  X{3}
  X{4}
  %}

% obtain circumcenter point
  A = [(X{1}-X{2})';
       (X{1}-X{3})';
       (X{1}-X{4})'];
  b = 0.5*[ X{1}'*X{1} - X{2}'*X{2};
            X{1}'*X{1} - X{3}'*X{3};
            X{1}'*X{1} - X{4}'*X{4}];

  tau1 = A\b;

% obtain shift  
  %tau = xi*tau1 + (1-xi)*tau0;

% shift  
  %R = R + tau;

  R = (1-xi)*tau0 + xi*tau1;

  
  
% echo bond length
%{
  fprintf(1,'bond length\n')
  for i = 1 : 4
  fprintf(1,'%d\t%.4f\n',i,sqrt((X{i}-R)'*(X{i}-R)));
  end
%}
  



end