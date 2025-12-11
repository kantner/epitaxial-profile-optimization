function [par] = set_J0_parameters(optimization_objective, par)  
    % set parameters for
    %      E1 + a1*M + a2*S  +  b1*nu + b2*sigma
    % J0 = -------------------------------------
    %      E2 + a3*M + a4*S  +  b3*nu + b4*sigma
    %
    % where M = <E_VS> and S = sqrt(Var(E_VS)) and E_1,2 are fixed
    % reference energies (for non-dimensionalization/ scaling)

    switch(optimization_objective)
      case 1
      % minimize ratio std(E_VS)/mean(E_VS)
        par.opt.E = [0 0]*par.units.meV; % E1, E2
        par.opt.a = [0 1 1 0]; % a1 (M), a2 (S), a3 (M), a4 (S)
        par.opt.b = [0 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)
      case 2
      % maximize mean(E_VS)
        par.opt.E = [1 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 1 0]; % a1 (M), a2 (S), a3 (M), a4 (S)
        par.opt.b = [0 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)
      case 3
      % minimize std(E_VS)
        par.opt.E = [0 1]*par.units.meV; % E1, E2
        par.opt.a = [0 1 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)
        par.opt.b = [0 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)
      case 4
      % maximize std(E_VS)
        par.opt.E = [1 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 1]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)
      case 5
      % minimize mean(E_VS)
        par.opt.E = [0 1]*par.units.meV; % E1, E2
        par.opt.a = [1 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)
      case 6
      % minimize mean(E_VS)/std(E_VS)
        par.opt.E = [0 0]*par.units.meV; % E1, E2
        par.opt.a = [1 0 0 1]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)        

      case 11
      % minimize mean(E_VS)/nu
        par.opt.E = [0 0]*par.units.meV; % E1, E2
        par.opt.a = [1 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 0 1 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)  
      case 12
      % minimize sigma/nu
        par.opt.E = [0 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 1 1 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)  
      case 13
      % minimize 1/nu
        par.opt.E = [1 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 0 1 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)    
      case 14
      % minimize sigma
        par.opt.E = [0 1]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 1 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)    
      case 15
      % minimize nu
        par.opt.E = [0 1]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [1 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)    
      case 16
      % minimize 1/sigma
        par.opt.E = [1 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [0 0 0 1]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)    
      case 17
      % minimize nu/sigma
        par.opt.E = [0 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 0 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [1 0 0 1]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)    
      case 18
      % minimize nu/mean
        par.opt.E = [0 0]*par.units.meV; % E1, E2
        par.opt.a = [0 0 1 0]; % a1 (M), a2 (S), a3 (M), a4 (S)  
        par.opt.b = [1 0 0 0]; % b1 (nu), b2 (sigma), b3 (nu), b4 (sigma)    

        
      otherwise
        error('invalid option')
    end

end