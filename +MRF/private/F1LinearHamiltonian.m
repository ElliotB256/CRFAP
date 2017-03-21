function [ H ] = F1LinearHamiltonian( t, Zs, RFs, gFuBB, phase )
%F1LINEARHAMILTONIAN Evaluates the Hamiltonian for F=1 atoms at the specified
%time. All MRF components are linearly polarised. There is currently no
%support for 3D potentials.
% 
% gFuBB: the amplitude of the rf field expressed in MHz

if nargin < 5
    phase = 0;
end

RFs = RFs / 2 / pi; % t is 2*pi*time

% Don't use input parser - this function must be called frequently. As a
% result, make it as lean as possible.
c = gFuBB;

H0 = [ 
       Zs,      0,      0,       ;
        0,      0,      0,       ;
        0,      0,    -Zs,       ;
     ];
         
% Assume: All rf components are circ pol. \alpha_i, \beta_i \to \alpha, \beta

% Calculate coherence terms
st    = sum( c .* sin(t .* 2 .* pi .* RFs + phase), 1) / (2.^0.5);
% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above. This is then divided by 2 to give 1/sqrt(2).

Hc = [     ...
         0,     st,      0,       ;
        st,      0,     st,       ;
         0,     st,      0,       ;
     ];
    
H = Hc + H0;

end

