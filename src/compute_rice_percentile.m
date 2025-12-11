function [x_alpha] = compute_rice_percentile(alpha, nu, sigma)

% cumulative pdf for integration to x_max
  cdf = @(x_max) quadgk(@(x) rician_pdf(x, nu, sigma), 0, x_max, 'AbsTol',1E-14); 

% define function to obtain upper bound such that cdf(x) = alpha  
  f = @(x) cdf(x) - alpha;
  
% compute mean of rice to initialize search
  [m, v] = rician_mean_and_variance(nu, sigma);

  x_alpha = fzero(@(x) f(x), m);

% plot
  DEBUG = 0;
  if DEBUG == 1
  figure(1342354); clf; hold all;

  x_range = linspace(0,m+3*sqrt(v));
  plot(x_range, rician_pdf(x_range, nu, sigma));
  
  xline(m,'b-','DisplayName','mean')


  xline(x_alpha,'r-','DisplayName','\alpha-percentile')
  legend()
  box on
  end

end