function [F] = nonlocal_well_F_0_s(X1, X2)
% function entering nonlocal pseudopotential for square s-well (l=0)
% nondimensionalized input X_i = K_i * R0 (where R0 is the well width)

% mean and diff
  dX = X2 - X1;
  Xm = (X1 + X2)/2;

% thresholds
  X_small  = 1E-12;

% case distinction
  % (1) if X1 is small but X2 is not
  % (2) if X2 is small but X1 is not
  % (3) if both X1 and X2 are nearly equal and small
  % (4) if both X1 and X2 are nearly equal but not small
  % (5) if both X1 + X2 is small but dX is not
  % (6) none of the above


  if abs(X1) < X_small && abs(X2) > X_small
    % case (1)
      eps = X1;
      X   = X2;

      Xsq   = X*X;
      sincX = sin(X)/X;
      cosX  = cos(X);

      F = (sincX - cosX)/(Xsq) ...
          + ((1 - 0.5*Xsq)*sincX + (1/6 * Xsq - 1)*cosX)/(Xsq*Xsq) * eps*eps;

  elseif abs(X2) < X_small && abs(X1) > X_small
    % case (2)
      eps = X2;
      X   = X1;

      Xsq   = X*X;
      sincX = sin(X)/X;
      cosX  = cos(X);

      F = (sincX - cosX)/(Xsq) ...
          + ((1 - 0.5*Xsq)*sincX + (1/6 * Xsq - 1)*cosX)/(Xsq*Xsq) * eps*eps;

  elseif abs(dX) < X_small && abs(Xm) < X_small
    % case (3)
      F = 1/3 - 1/30 * (X1*X1 + X2*X2);

  elseif abs(dX) < X_small && abs(Xm) >= X_small      
    % case 4
      Xm_sq  = Xm*Xm;
      cosXm  = cos(Xm);
      sincXm = sin(Xm)/Xm;

      F = 1/(2*Xm_sq) * (1 - sincXm * cosXm) ...
          + 1/(8 * Xm_sq*Xm_sq) * (1 - 2/3*Xm_sq - cosXm * sincXm) * dX*dX;

  elseif abs(Xm) < X_small && abs(dX) >= X_small
    % case 5
      dX_sq  = dX*dX;
      sincdX = sin(dX)/dX;

      F = 2/(dX_sq) * (1 - sincdX) ...
          + 8/(dX_sq*dX_sq) * (1 - 1/6*dX_sq - sincdX) * Xm*Xm;
  else
    % case 6
      F = (sin(X1)/X1 * cos(X2) - sin(X2)/X2 * cos(X1))/( (X1-X2)*(X1+X2) );
  end





end