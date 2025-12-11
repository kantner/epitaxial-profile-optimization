function [x_ww] = x_wiggle_well(x_ampl, q, phi, par)

  %x_ww = x_ampl * 0.5*(1 + cos(q*par.z + phi)) .* par.QW_indicator;

  mode = 2;
  switch(mode)
    case 1
      trig_func = @(x) sin(x);
    case 2
      trig_func = @(x) cos(x);
    case 3
      trig_func = @(x) -sin(x);
    case 4
      trig_func = @(x) -cos(x);
  end


  x_ww = x_ampl * 0.5*(1 + trig_func(q*par.z + phi)) .* par.QW_indicator;
end