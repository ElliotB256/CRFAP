function [ H ] = F2Hamiltonian( t, Zs, RFs, gFuBB, theta )
%F2HAMILTONIAN Evaluates the Hamiltonian for F=2 atoms at the specified
%time. All MRF components are circularly polarised.

if nargin < 5
    theta = 0;
end

RFs = RFs / 2 / pi; % t is 2*pi*time

% Don't use input parser - this function must be called frequently. As a
% result, make it as lean as possible.
e = exp(1);
c = gFuBB;

H0 = [  2 * Zs,     0,      0,      0,      0       ;
             0,    Zs,      0,      0,      0       ;
             0,     0,      0,      0,      0       ;
             0,     0,      0,    -Zs,      0       ;
             0,     0,      0,      0,  -2*Zs       ];
         
% Assume: All rf components are circ pol. \alpha_i, \beta_i \to \alpha, \beta
alpha = (cos(theta) - 1)/2;
beta  = (cos(theta) + 1)/2;

% Calculate coherence terms
cp    = sum( c .* ( alpha * e.^( 1i .* t .* 2 .* pi .* RFs ) + beta * e.^(-1i .* t .* 2 .* pi .* RFs ) ), 1) / 2;
cn    = sum( c .* ( alpha * e.^(-1i .* t .* 2 .* pi .* RFs ) + beta * e.^( 1i .* t .* 2 .* pi .* RFs ) ), 1) / 2;

Hc = [     ...
           0,           2*cp,           0,              0,      0       ;
        2*cn,              0,   6.^0.5*cp,              0,      0       ;
           0,      6.^0.5*cn,           0,      6.^0.5*cp,      0       ;
           0,              0,   6.^0.5*cn,              0,   2*cp       ;
           0,              0,           0,           2*cn,      0       ;
        ];
    
% oscillating component of rf field parallel to local field
% pc = parallel component
pc = sum(c .* sin(theta) .* cos(2 .* pi .* RFs .* t), 1);
Hp = [ ...
    2 .* pc,        0,      0,      0,      0       ;
          0,       pc,      0,      0,      0       ;
          0,        0,      0,      0,      0       ;
          0,        0,      0,    -pc,      0       ;
          0,        0,      0,      0,  -2*pc       ;
          ];
    
H = Hc + Hp + H0;

end

