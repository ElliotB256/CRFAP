function [ w ] = getTrapFreq( x, pot, massAmu )
%GETTRAPFREQ Gets trap frequencies from the given potential.
% x: positions in um
% pot: potential in MHz
% massAmu: mass of species in amu

import Constants.*

%convert the trap fit to realistic frequencies.
factor = ((pi ) *1e6 / (massAmu) * hbar/amu).^0.5 * 1e6 /(pi);

k = performHarmonicFit(x, pot);
w = k.^0.5 * factor;


end

