function [ delta ] = MRFPropTimeEv( t, ds, RF, BRF, B )
%TIMEEV 

% coherence term
c = sum(BRF .* cos(RF .* t), 1) / 2.^0.5;

% calculate product of Hamiltonian and ds in column form
% This looks ugly!

H = [ B, c, 0; c, 0, c ; 0, c, -B ];
i = [ 1 2 3; 4 5 6; 7 8 9;]';

ds = ds(i);
delta = -1i * H * ds;
delta = delta(:);

end

