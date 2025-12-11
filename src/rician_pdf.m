function f = rician_pdf(x, nu, sigma)

  sigma_sq = sigma*sigma;
  f = x/sigma_sq .* exp( -0.5*(x.^2 + nu.^2)./sigma_sq) .* besseli(0, x*nu/sigma_sq);

end