function [ H ] = Hamiltonian( B0, RFs, Rabi )
%HAMILTONIAN Calculates the Hamiltonian for a MRF system described by the
%specified parameters.
% Syntax: Hamiltonian( B0, RFs, Rabi )
%  Bs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)

c = Rabi*2;

H = @(t) [ B0 , sum(c .* cos(RFs .* t), 1) / 2.^0.5, 0;
    sum(c .* cos(RFs .* t), 1) / 2.^0.5, 0, sum(c .* cos(RFs .* t), 1) / 2.^0.5;
    0, sum(c .* cos(RFs .* t), 1) / 2.^0.5, -B0 ];

end

