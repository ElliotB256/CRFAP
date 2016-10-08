function [ delta ] = PropTimeEv( t, ds, H )
%PropTimeEv Calculates the differential at specified time, required for
%propagator calculation. H is a function handle that accepts one argument,
%current time, and returns a 3x3 matrix.

h = H(t);
i = [ 1 2 3; 4 5 6; 7 8 9;]';
ds = ds(i);
delta = -1i * h * ds;
delta = delta(:);

end