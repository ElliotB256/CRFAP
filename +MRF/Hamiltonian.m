function [ H ] = Hamiltonian( Zs, RFs, Rabi, varargin )
%HAMILTONIAN Calculates the Hamiltonian for a MRF system described by the
%specified parameters.
% Syntax: Hamiltonian( Zs, RFs, Rabi )
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)
%  F: Specify F=1 or F=2 system.
%  theta: used to calculate circ polarised AP at arbitrary position.

p = inputParser;
addRequired(p,'Zs',@isnumeric);
addRequired(p,'RFs',@isnumeric);
addRequired(p,'Rabi',@isnumeric);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));
addParameter(p, 'theta', 0, @(x) x >= 0 && x <= pi);

parse(p,Zs, RFs, Rabi, varargin{:});

theta = p.Results.theta;

switch p.Results.F
    case 1
        % Circ pol with arb rotation:
        H = @(t) F1Hamiltonian(t, Zs, RFs, Rabi, theta);
    case 2
        % Circ pol with arb rotation:
        H = @(t) F2Hamiltonian(t, Zs, RFs, Rabi, theta);
    otherwise
        error('Unsupported value of F specified.');
end

end

