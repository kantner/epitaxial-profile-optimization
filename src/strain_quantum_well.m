function [eps_QW] = strain_quantum_well(x_subs, par)
% compute strain tensor of biaxially strain Si on Si_{1-x_subs} Ge_{x_subs}
    a0_subs = par.a0_SiGe(x_subs);

    eps_QW      = zeros(3,3);  

    eps_QW(1,1) = (a0_subs - par.a0)/par.a0;
    eps_QW(2,2) = eps_QW(1,1);
    eps_QW(3,3) = -par.C12/par.C11 * (eps_QW(1,1) + eps_QW(2,2));

end