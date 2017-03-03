function [ w ] = getTrapFreq( x, pot, massAmu )
%GETTRAPFREQ Gets trap frequencies from the given potential.
% x: positions in um
% pot: potential in MHz
% massAmu: mass of species in amu
% Syntax: getTrapFreq( x, pot, massAmu)
% Output is in Hz

import Constants.*
import Util.*

%convert the trap fit to realistic frequencies.
factor = ((pi ) *1e6 / (massAmu) * hbar/amu).^0.5 * 1e6 /(pi);

% U = m omega^2 z^2 / 2
% ((MHz) * 1e6) * h = (massAmu * amu) * ( 2 pi * f ) ^2 * (z in um * 1e6)^2

% Checked factor 29/11/2016 - correct.

k = performHarmonicFit(x, pot);
w = k.^0.5 * factor;

end

