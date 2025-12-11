function [r] = two_loop_recursion(grad_f, s, y)
 % Implementation of two-loop recursion for computation of the LBFGS search direction by evaluating the action of inverse quasi-Hessian on the gradient
 % INPUT
 %   s_k ... history of updates          s_k = x_k+1 - x_k = alpha*p_k
 %   y_k ... history of gradient updates y_k = gradf_k+1 - gradf_k
 % OUTPUT
 %   r = invH*gradf (ascent direction!)
   
  % system size
    %n = length(grad_f);
    n = size(s, 1);

  % memory depth
    m = size(s, 2);
    
  % allocate memory
    alpha = zeros(1,m);
    rho   = zeros(1,m);

  % set initial inverse Hessian (only for m>0)    
    if m > 0
      set_scaling = 1;
      switch(set_scaling)
        case 0 % no scaling, use invH0 = Id
          %invH0 = speye(n);
          invH0 = 1;% unit matrix
        case 1 % set scaled invH0 matrix according to Nocedal-Wright formula (6.20)     
          %invH0 = (y(:,1)' * s(:,1))/(y(:,1)'*y(:,1)) * speye(n);
          invH0 = (y(:,m)' * s(:,m))/(y(:,m)'*y(:,m));
        otherwise
          error('not implemented')
      end

    else
      invH0 = 1; % unit matrix
    end
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % start two-loop recursion
    q = grad_f;

  % right product
    for i = m : -1 : 1
      rho(i)   = 1/(y(:,i)' * s(:,i));
      alpha(i) = rho(i) * (s(:,i)' * q);
      q        = q - alpha(i) * y(:,i);
    end

  % center product  
    r = invH0 * q;

  % left product  
    for i = 1 : 1 : m
      beta = rho(i) * (y(:,i)' * r);
      r    = r + s(:,i)*(alpha(i)-beta);
    end

  %%%%%%%%%%%%%%%%%  
  % check
  %{
    % curvature conditions
      curv = zeros(m,1);
      for i = 1 : m
        curv(i) = y(:,i)' * s(:,i);
      end
      curv
  %}


    if sum(isnan(r)) > 0      
      error('r contains NaN')
    end

    if sum(isinf(r)) > 0
      error('r contains Inf')
    end

end