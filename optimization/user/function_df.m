function [df, out] = function_df(x, in, par)
  %{
model = par.model;

  switch(model)
    case 1
      df = df_rosenbrock(x, par);
    case 2      
      df = df_styblinksi_tang(x, par);
    otherwise
      error('not implemented')
  end
  %}
  [out] = cost_gradient(x, in, par);
  df = out.DJ;

  if sum(isnan(df))
    error('df contains NaN')
  end

  if sum(isinf(df))
    error('df contains Inf')
  end

end