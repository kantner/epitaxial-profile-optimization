function [] = optimization_map_scan_eval(results, par)

% obtain size
  nx = length(results.x_mean_range);
  nk = length(results.kc_range);


  if min(nk,nx) == 1 
    optimization_map_scan_eval_1D(results, par)    
  else
    optimization_map_scan_eval_2D(results, par)
  end


end