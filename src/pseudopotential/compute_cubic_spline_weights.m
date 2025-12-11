function [w] = compute_cubic_spline_weights(par)
% Compute cubic spline interpolation parameters (w1 ... w16)
% V_loc (q) = w1*q^3 + w2*q^2 + w3*q + w4
% in sections (i)         0  < q < sqrt(3)     (parameters w1...w4)
%             (ii)  sqrt(3)  < q < sqrt(8)     (parameters w5...w8)
%             (iii) sqrt(8)  < q < sqrt(11)    (parameters w9...w12)
%             (iv)  sqrt(11) < q < 3 k_F       (parameters w13...w16)
% BCs:        V_loc(0)   = -2/3*EF
%             V_loc(3 k_F)  = 0
%             V_loc'(3 k_F) = 0
%      (v1)   V_loc'(0) = 0
%      (v2)   V_loc''(0)  = 0
%
% Here q and k_F is assumed in units of 2*pi/a0 (nondimensionalized).
%
% There are two variants of the BC setting either V_loc'(0)=0 or V_loc''(0)=0
  BC_mode = par.pp.BC_mode;

  if ismember(BC_mode, [1 2]) ~= 1
    error('BC_mode must be 1 or 2')
  end

% import
  V_0      = par.pp.V_0      * 1/par.units.Ry;
  V_sqrt3  = par.pp.V_sqrt3  * 1/par.units.Ry;
  V_sqrt8  = par.pp.V_sqrt8  * 1/par.units.Ry;
  V_sqrt11 = par.pp.V_sqrt11 * 1/par.units.Ry;

  kF       = par.pp.kF / (2*pi/par.a0);

% create matrix and rhs vector
  A = zeros(16, 16);
  b = zeros(16, 1);


 %% match local potential values (3 conditions)

  % condition
    idx = 1;
    q   = sqrt(3);
    
    A(idx,1) = q^3;
    A(idx,2) = q^2;
    A(idx,3) = q^1;
    A(idx,4) = 1;
    b(idx)   = V_sqrt3;
    
  % condition
    idx = 2;
    q   = sqrt(8);
    
    A(idx,5) = q^3;
    A(idx,6) = q^2;
    A(idx,7) = q^1;
    A(idx,8) = 1;
    b(idx)   = V_sqrt8;
    
  % condition
    idx = 3;
    q   = sqrt(11);
    
    A(idx,9)  = q^3;
    A(idx,10) = q^2;
    A(idx,11) = q^1;
    A(idx,12) = 1;
    b(idx)    = V_sqrt11;

 %% internal matching at sqrt(3)
    q   = sqrt(3);
  
  % continuity
    idx = 4;
    
    A(idx,1) = q^3;
    A(idx,2) = q^2;
    A(idx,3) = q^1;
    A(idx,4) = 1;
    
    A(idx,5) = -q^3;
    A(idx,6) = -q^2;
    A(idx,7) = -q^1;
    A(idx,8) = -1;
    
  % differentiabilty
    idx = 5;
    
    A(idx,1) = 3*q^2;
    A(idx,2) = 2*q^1;
    A(idx,3) = 1;
    
    A(idx,5) = -3*q^2;
    A(idx,6) = -2*q^1;
    A(idx,7) = -1;
  
  % curvature
    idx = 6;
    
    A(idx,1) = 2*3*q^1;
    A(idx,2) = 1*2;
    
    A(idx,5) = -2*3*q^1;
    A(idx,6) = -1*2;

 %% internal matching at sqrt(8)
    q   = sqrt(8);
    
  % continuity
    idx = 7;
    
    A(idx,5) = q^3;
    A(idx,6) = q^2;
    A(idx,7) = q^1;
    A(idx,8) = 1;
    
    A(idx,9)  = -q^3;
    A(idx,10) = -q^2;
    A(idx,11) = -q^1;
    A(idx,12) = -1;
       
  % differentiabilty
    idx = 8;
    
    A(idx,5) = 3*q^2;
    A(idx,6) = 2*q^1;
    A(idx,7) = 1;
    
    A(idx,9) = -3*q^2;
    A(idx,10) = -2*q^1;
    A(idx,11) = -1;
    
  % curvature
    idx = 9;
    
    A(idx,5) = 2*3*q^1;
    A(idx,6) = 1*2;
    
    A(idx,9) = -2*3*q^1;
    A(idx,10) = -1*2;

 %% internal matching at sqrt(11)
    q   = sqrt(11);

  % continuity
    idx = 10;
    
    A(idx,9)  = q^3;
    A(idx,10) = q^2;
    A(idx,11) = q^1;
    A(idx,12) = 1;
    
    A(idx,13) = -q^3;
    A(idx,14) = -q^2;
    A(idx,15) = -q^1;
    A(idx,16) = -1;

  % differentiabilty
    idx = 11;
    
    A(idx,9)  = 3*q^2;
    A(idx,10) = 2*q^1;
    A(idx,11) = 1;
    
    A(idx,13) = -3*q^2;
    A(idx,14) = -2*q^1;
    A(idx,15) = -1;

  % curvature
    idx = 12;
    
    A(idx,9)  = 2*3*q^1;
    A(idx,10) = 1*2;
    
    A(idx,13) = -2*3*q^1;
    A(idx,14) = -1*2;

 %% BC for V at q=0
    q   = 0;

  % V(0) = V_0
    idx = 13;
    
    A(idx,1) = q^3;
    A(idx,2) = q^2;
    A(idx,3) = q^1;
    A(idx,4) = 1;
    
    b(idx) = V_0;

 %% BC for 1st or 2nd derivative of V at q=0
    idx = 14;

    switch(BC_mode)
      case 1 % V'(0) = 0
        A(idx,1) = 3*q^2;
        A(idx,2) = 2*q^1;
        A(idx,3) = 1;
        b(idx)   = 0;
      case 2 % V''(0) = 0
        A(idx,1) = 2*3*q;
        A(idx,2) = 1*2;
        b(idx)   = 0;
      otherwise
        error('BC_mode must be 1 or 2')
    end

 %% BCs at q=3 k_F
    q   = 3*kF;

  % V(3 k_F) = 0
    idx = 15;
    
    A(idx,13) = q^3;
    A(idx,14) = q^2;
    A(idx,15) = q^1;
    A(idx,16) = 1;
    
    b(idx) = 0;

  % V'(3k_F) = 0
    idx = 16;
    
    A(idx,13) = 3*q^2;
    A(idx,14) = 2*q^1;
    A(idx,15) = 1;
    
    b(idx) = 0;


%%%%%%%%%%%%%%%%%%%%
 %% solve
    w = A\b;

end








