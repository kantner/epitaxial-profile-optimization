function [phi] = linesearch_phi(x, alpha, p, par)
% evaluate line search function
% phi(alpha) = J(x + alpha * p)




x = x + alpha * p;

compute_derivatives = 1;
out = cost_functional(x, compute_derivatives, par);

phi = out.J_tot;

end