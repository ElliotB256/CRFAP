function [ eigF2 ] = GetQuasiEnergies( Bs, RF, BRF )
%GETQUASIENERGIES Calculates the quasi-energies for the given multi-RF
%field over the specified range of zeeman energy splittings.

% Note: At present, the Bs are specified in MHz. This means the factor gF
% uB has been incoporated into the B specified here. Think of it more as
% the energy splitting of the bare states.

eigF2 = zeros(3, length(Bs));
periodicity = 2*pi/MRF.GetFundamental(RF);
for i=1:length(Bs)
    B = Bs(i);
    H = MRF.Hamiltonian(B, RF, BRF);
    U = MRF.Propagator(H, periodicity);
    eigF2(:,i) = sort(angle(eig(U)))/periodicity;
end

end

