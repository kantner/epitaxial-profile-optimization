function [F] = generate_filter(K, K_c, filter_type)
% generate frequency-domain filter function
  
  switch(filter_type)
    case 0 % no filtering
      F = ones(size(K));
    case 1 % hard cutoff
      F = ones(size(K));
      F(abs(K)>K_c) = 0;
    case 2 % exponential filter
      F = exp(-log(2) * abs(K)/K_c);
    case 3 % gaussian filter
      F = exp(-log(2) * (K/K_c).^2 );      
    case 4 % bessel filter
      F = bessel_filter(85, K, K_c);
    case 5 % sigmoid (fermi)
      dK = 0.1E9;
      F  = 1./(exp((abs(K)-K_c)/dK) + 1);
    case 6 % sigmoid (erf)
      dK = 0.1E9;
      X  = (abs(K)-K_c)/dK;
      F  = 0.5*(1+erf(-X));
      
    otherwise
      error('not implemented');
  end



% plot
%{
  if filter_type ~= 0

    figure(99999); clf; hold all;
      [~,srt] = sort(K);
      plot(K(srt),F(srt),'k-')
      yline(0.5,'r--')
      xline(+K_c,'b--')
      xline(-K_c,'b--')
      xlabel('k')
      ylabel('F(k)')
      box on
      title('filter funtion')
      
  end
%}

end