function [ mu ] = chemicalPotential( omega, N, a, m )
%CHEMICALPOTENTIAL Calculates the chemical potential of atoms in the
%specified trap. The dimensionality used depends on the dimensionality of
%omega. Formulae from Merloti, N. J. Phys, 2012.
%
% Omega             : trap frequencies, 2pi x MHz
% 
% N                 : Number of atoms
%
% a                 : scattering length, Bohr radii
%
% m                 : atomic mass, amu
%
% Returns: chemical potential in MHz

import Constants.*;

if length(omega) ~= 2 && length(omega) ~= 3
    error('Only 2D and 3D cases are supported');
end

omegaBar = geomean(omega);

% oscillator Lengths in Bohr radii

switch (length(omega))
    case 3        
        oscillatorLength = ( (hbar/amu) / ( m * omegaBar) ).^0.5 / bohr;
        mu = (omegaBar / (2 * pi)) / 2 * (15 * N * a / oscillatorLength) .^ (2/5);
    case 2
        oscillatorLength = ( (hbar/amu) / ( m * max(omegaBar(:))) ).^0.5 / bohr;
        mu = (2 * min(omegaBar(:)) / (2 * pi)) * (N * a / (2 * pi).^0.5 / oscillatorLength) .^ 1/2;
end

end

