function [ H ] = Hamiltonian( Zs, RFs, gFuBB, varargin )
%HAMILTONIAN Calculates the Hamiltonian for a MRF system described by the
%specified parameters.
% Syntax: Hamiltonian( Zs, RFs, Rabi )
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  gFuBB: vector of dressing RF Rabi frequencies (MHz)
%  F: Specify F=1 or F=2 system.
%  theta: used to calculate circ polarised AP at arbitrary position.

p = inputParser;
addRequired(p,'Zs',@(x) isnumeric(x) && size(x,2) == 1);
addRequired(p,'RFs',@(x) isnumeric(x) && size(x, 2) == 1);
addRequired(p,'Rabi',@(x) isnumeric(x) && size(x, 2) == 1);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));
addParameter(p, 'theta', 0, @(x) all(x >= 0 & x <= pi) && size(x, 2) == 1);

parse(p,Zs, RFs, gFuBB, varargin{:});

theta = p.Results.theta;

switch p.Results.F
    case 1
        % Circ pol with arb rotation:
        H = @(t) F1Hamiltonian(t, Zs, RFs, gFuBB, theta);
    case 2
        % Circ pol with arb rotation:
        H = @(t) F2Hamiltonian(t, Zs, RFs, gFuBB, theta);
    otherwise
        error('Unsupported value of F specified.');
end

end

