function [x] = apply_filter(x, filter)
  x = ifft( filter .* fft(x));
end
