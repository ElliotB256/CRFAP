function [ epsilon ] = calculateEpsilon( zsf, BGrad, M )
%CALCULATEEPSILON Calculates the value of epsilon, which is a ratio of
%gravity to magnetic confinement strength defined in Merloti, N. J. Phys,
%2012.
%
% epsilon = M g / (2 m_F hbar alpha)
%
% zsf   : Zeeman split factor, g_f \mu_B / hbar, in units of MHz/Gauss.
% 
% BGrad : Quadrupole gradient, Gauss/cm
% 
% M     : atomic mass
% 
% Returns: dimensionless constant

% get gpe for 1 micron to get gp in units of MHz/micron
mgh = (gpe(1, M));

% calculate magnetic energy, convert BGrad to microns
magneticEnergy = zsf * (BGrad * 1e-4);

epsilon = mgh / magneticEnergy / 2;

end

