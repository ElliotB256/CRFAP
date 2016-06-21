function [ h ] = calcH( B0, RFs, BRFs, t )
%CALCH Calculates Hamiltonian for given frequencies and time

p = sum(BRFs .* cos(RFs .* t), 1) / 2.^0.5;
h = [ B0, p, 0; p, 0, p; 0, p, -B0 ];

end