function [ pot ] = gpe( z, mass )
%GPE Gravitational potential energy in MHz
% mass: mass of atom in amu

amu = 1.661e-27; %kg
g = 9.81;

%convert g to microns
g = g * 1e-6;

% and now calculate gpe in MHz
h = 6.63e-34;
pot = z * mass * g * (amu / h) * 1e-6; %1e-6 convert from Hz to MHz

end

