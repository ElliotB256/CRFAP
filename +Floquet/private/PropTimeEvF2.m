function [ delta ] = PropTimeEvF2( t, ds, H )
%PropTimeEv Calculates the differential at specified time, required for
%propagator calculation. H is a function handle that accepts one argument,
%current time, and returns a 5x5 matrix. This version operates in a 5x5
%space used for F=2 evaluation.

h = H(t);
i = [ 1:5; 6:10; 11:15; 16:20; 21:25]';
ds = ds(i);
delta = -1i * h * ds;
delta = delta(:);

end