function [ H ] = Hamiltonian( B0, RFs, Rabi )
%HAMILTONIAN Returns the Hamiltonian for an F=1 system irradiated by the
%given RFs

c = Rabi*2;

%p = @(t) sum(BRFs .* cos(RFs .* t), 1) / 2.^0.5;
H = @(t) [ B0 , sum(c .* cos(RFs .* t), 1) / 2.^0.5, 0;
    sum(c .* cos(RFs .* t), 1) / 2.^0.5, 0, sum(c .* cos(RFs .* t), 1) / 2.^0.5;
    0, sum(c .* cos(RFs .* t), 1) / 2.^0.5, -B0 ];

end

