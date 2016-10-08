function [ H ] = Hamiltonian( Zs, RFs, Rabi, varargin )
%HAMILTONIAN Calculates the Hamiltonian for a MRF system described by the
%specified parameters.
% Syntax: Hamiltonian( Zs, RFs, Rabi )
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)
%  Optional parameters: the values of alpha, Beta, zeta / Rabi

p = inputParser;
addRequired(p,'Zs',@isnumeric);
addRequired(p,'RFs',@isnumeric);
addRequired(p,'Rabi',@isnumeric);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));

parse(p,Zs, RFs, Rabi, varargin{:});


% c = Rabi*2;
% H = @(t) [ B0 , sum(c .* cos(RFs .* t), 1) / 2.^0.5, 0;
%     sum(c .* cos(RFs .* t), 1) / 2.^0.5, 0, sum(c .* cos(RFs .* t), 1) / 2.^0.5;
%     0, sum(c .* cos(RFs .* t), 1) / 2.^0.5, -B0 ];

% change Hamiltonian to use circ pol. [0, exp(iwt), 0; exp(-iwt), 0, exp(iwt); 0, exp(iwt), 0 ];
e = exp(1);
c = Rabi;

switch F
    case 1
        
        H = @(t) [ ...
            Zs, sum( c .* e .^ (-RFs .* 1i .* t), 1) / 2.^0.5, 0;
            sum( c .* e .^ (RFs .* 1i .* t), 1) / 2.^0.5, 0, sum( c .* e .^ (-RFs .* 1i .* t), 1) / 2.^0.5;
            0, sum( c .* e .^ (RFs .* 1i .* t), 1) / 2.^0.5, -Zs;
            ];
        
    case 2
        
        % s = sign of exponent
        omega = @(t,s) sum( c .* e .^ (s*RFs .* 1i .* t), 1) / 2.^0.5;
        H = @(t) [ ...
                   2*Zs,    omega(t,-1),              0,              0,              0;
            omega(t, 1),             Zs,    omega(t,-1),              0,              0;
                      0,     omega(t,1),              0,    omega(t,-1),              0;
                      0,              0,     omega(t,1),            -Zs,    omega(t,-1);
                      0,              0,              0,     omega(t,1),          -2*Zs;
                 ];
        
    otherwise
        error('Unsupported value of F specified.');
end

end

