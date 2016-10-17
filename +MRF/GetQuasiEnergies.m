function [ eigF2 ] = GetQuasiEnergies( Zs, RF, Rabi, varargin )
%GETQUASIENERGIES Calculates quasi-energies for the given multi-RF
%field over the specified range of zeeman energy splittings.
% Syntax: GetQuasiEnergies( Zs, RF, Rabi)
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  RF: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)

p = inputParser;
addRequired(p,'Zs',@isnumeric);
addRequired(p,'RF',@isnumeric);
addRequired(p,'Rabi',@isnumeric);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));
addParameter(p, 'theta', 0, @(x) x >= 0 && x <= pi);
parse(p,Zs, RF, Rabi, varargin{:});

switch p.Results.F
    case 1
        hs = 3;
    case 2
        hs = 5;
end

eigF2 = zeros(hs, length(Zs));
periodicity = 2*pi/MRF.GetFundamental(RF);

%periodicity = 1/MRF.GetFundamental(RF);
%warning('blah');

for i=1:length(Zs)
    B = Zs(i);
    H = MRF.Hamiltonian(B, RF, Rabi, 'F', p.Results.F, 'theta', p.Results.theta);
    U = MRF.Propagator(H, periodicity, p.Results.F);
    eigF2(:,i) = sort(angle(eig(U)))/periodicity;
end

end