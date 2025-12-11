function L = Laplacian(N)
% generate 1D finite difference laplacian

  periodic = 1;

  if periodic == 1
    L = -2*diag(sparse(ones(N,1))) ...
        + circshift(diag(sparse(ones(N,1))),+1) ...
        + circshift(diag(sparse(ones(N,1))),-1);

  else
    L = -2 * diag(sparse(ones(N,1))) ...
        + diag(sparse(ones(N-1,1)),+1) ...
        + diag(sparse(ones(N-1,1)),-1);
  end

end