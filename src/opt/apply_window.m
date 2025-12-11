function [x] = apply_window(x, window)
  x = window .* x;
end
