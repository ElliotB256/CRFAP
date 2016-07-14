function [ pot ] = gpe( z, mass )
%GPE Gravitational potential energy in MHz
% mass: mass of atom in amu
% z: microns
% SYNTAX: gpe( z, mass )

%convert z to metres
z = z * 1e-6;

% and now calculate gpe in MHz
pot = z * mass * Constants.g * (Constants.amu / Constants.h) * 1e-6; %1e-6 convert from Hz to MHz

end