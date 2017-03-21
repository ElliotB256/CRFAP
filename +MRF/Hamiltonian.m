function [ H ] = Hamiltonian( Zs, RFs, gFuBB, varargin )
%HAMILTONIAN Calculates the Hamiltonian for a MRF system described by the
%specified parameters.
% Syntax: Hamiltonian( Zs, RFs, Rabi )
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  phase: phase vector of the dressing rf components
%  gFuBB: vector of dressing RF Rabi frequencies (MHz)
%  F: Specify F=1 or F=2 system.
%  theta: used to calculate circ polarised AP at arbitrary position.
%  polarisation: specify whether to use 'circ' or 'lin' polarisation.

p = inputParser;
addRequired(p,'Zs',@(x) isnumeric(x) && size(x,2) == 1);
addRequired(p,'RFs',@(x) isnumeric(x) && size(x, 2) == 1);
addRequired(p,'Rabi',@(x) isnumeric(x) && size(x, 2) == 1);
addParameter(p, 'phase', 0, @(x) all(x >= 0 & x <= 2*pi) && size(x, 2) == 1);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));
addParameter(p, 'theta', 0, @(x) all(x >= 0 & x <= pi) && size(x, 2) == 1);
addParameter(p, 'polarisation', 'circ', @(x) ismember(x, {'circ', 'lin'}));


parse(p,Zs, RFs, gFuBB, varargin{:});

theta = p.Results.theta;

switch p.Results.polarisation
    case 'circ'
        
        switch p.Results.F
            case 1
                % Circ pol with arb rotation:
                H = @(t) F1Hamiltonian(t, Zs, RFs, gFuBB, theta, p.Results.phase);
            case 2
                % Circ pol with arb rotation:
                H = @(t) F2Hamiltonian(t, Zs, RFs, gFuBB, theta, p.Results.phase);
            otherwise
                error('Unsupported value of F specified.');
                
        end
    case 'lin'
        
        if p.Results.F ~= 1
            error('Unsupported value of F specified for linear polarisation');
        end
        
        
end

end

