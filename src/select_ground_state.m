function idx_gnd = select_ground_state(S, E, psi_ref, x, par)
% selection of ground state in potential with electric field
% heuristic takes overlap with QW and energetic distance to confinement
% potential minimum into account

  %%%%%%%%%%%%%%%%%%%%
  % compute overlap with reference function  
    overlap = zeros(par.neigs,1);
    for i = 1 : par.neigs
      overlap(i) = abs(psi_ref' * S(:,i))/(psi_ref'*psi_ref);
    end       
    
  %%%%%%%%%%%%%%%%%%%%
  % obtain distance to reference energy

  % total potential energy
    U_X = par.dEc * x;
    U   = par.U_QW + par.U_F + U_X;

  % obtain minimum energy in QW domain and set as reference
    threshold = 0.95;
    %idx       = par.QW_indicator > threshold;
    idx = find(par.QW_indicator > threshold * max(par.QW_indicator));
    %if isempty(idx)
    %  U_min = 0;
    %else
      U_min = min(U(idx));
    %end
    %U_min

  % compute energy distance
    energy_dist = zeros(par.neigs,1);
    for i = 1 : par.neigs
      energy_dist(i) = abs(E(i)-U_min)/par.energy_scale;
    end

  %%%%%%%%%%%%%%%%%%%%%
  % selection criterion  
    %criterion = overlap./sqrt(energy_dist);
    criterion = overlap./(energy_dist.^1.0);

    %[overlap, energy_dist, criterion]
          
  % find maximum of criterion  
    idx_gnd = find(criterion == max(criterion));
    
end
    