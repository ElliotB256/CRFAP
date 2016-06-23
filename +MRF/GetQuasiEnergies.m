function [ eigF2 ] = GetQuasiEnergies( Bs, RF, Rabi )
%GETQUASIENERGIES Calculates quasi-energies for the given multi-RF
%field over the specified range of zeeman energy splittings.
% Syntax: GetQuasiEnergies( Bs, RF, Rabi)
%  Bs: energy splitting of the undressed Zeeman states in MHz.
%  RF: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)

eigF2 = zeros(3, length(Bs));
periodicity = 2*pi/MRF.GetFundamental(RF);
for i=1:length(Bs)
    B = Bs(i);
    H = MRF.Hamiltonian(B, RF, Rabi);
    U = MRF.Propagator(H, periodicity);
    eigF2(:,i) = sort(angle(eig(U)))/periodicity;
end

end