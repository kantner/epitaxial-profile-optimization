%% some physical constants
 % all in SI basis unit

   planckConstant        = 6.62607015 * 1E-34;    % (h) exact                       
   reducedPlanckConstant = planckConstant/(2*pi); % (hbar) exact
   elementaryCharge      = 1.602176634 * 1E-19;   % exact
   boltzmannConstant     = 1.380649 * 1E-23;      % exact                        
   vacuumSpeedOfLight    = 299792458;             % exact
  
   vacuumDielectricConstant = 8.8541878128 * 1E-12;  
   magneticConstant      = 1/(vacuumSpeedOfLight^2*vacuumDielectricConstant);
                             
   electronMass          = 9.1093837015 * 1E-31;
                       
   bohrRadius            = 4*pi*vacuumDielectricConstant*reducedPlanckConstant^2/(electronMass*elementaryCharge^2);
   bohrMagneton          = elementaryCharge*reducedPlanckConstant/(2*electronMass);

   %rydberg               = elementaryCharge^2/(8*pi*vacuumDielectricConstant*bohrRadius);
   rydberg               = 1/(2*electronMass) * (reducedPlanckConstant/bohrRadius)^2;
   hartree               = 2*rydberg;