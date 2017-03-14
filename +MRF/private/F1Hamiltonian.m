function [ H ] = F1Hamiltonian( t, Zs, RFs, gFuBB, theta, phase )
%F1HAMILTONIAN Evaluates the Hamiltonian for F=1 atoms at the specified
%time. All MRF components are circularly polarised.

if nargin < 5
    theta = 0;
end

if nargin < 6
    phase = 0;
end

RFs = RFs / 2 / pi; % t is 2*pi*time

% Don't use input parser - this function must be called frequently. As a
% result, make it as lean as possible.
e = exp(1);
c = gFuBB;

H0 = [ 
       Zs,      0,      0,       ;
        0,      0,      0,       ;
        0,      0,    -Zs,       ;
     ];
         
% Assume: All rf components are circ pol. \alpha_i, \beta_i \to \alpha, \beta
alpha = (cos(theta) - 1)/2;
beta  = (cos(theta) + 1)/2;

% Calculate coherence terms
cp    = (2.^0.5) * sum( c .* ( alpha .* e.^(1i .* (t .* 2 .* pi .* RFs + phase)) + beta .* e.^(-1i .* (t .* 2 .* pi .* RFs + phase)) ), 1) / 2;
cn    = (2.^0.5) * sum( c .* ( alpha .* e.^(-1i .* (t .* 2 .* pi .* RFs + phase)) + beta .* e.^(1i .* (t .* 2 .* pi .* RFs + phase)) ), 1) / 2;

% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above.
Hc = [     ...
         0,     cp,      0,       ;
        cn,      0,     cp,       ;
         0,     cn,      0,       ;
     ];
    
% oscillating component of rf field parallel to local field
% pc = parallel component
pc = sum(c .* sin(theta) .* cos(2 .* pi .* RFs .* t + phase), 1);
Hp = [ ...
          pc,      0,      0,       ;
           0,      0,      0,       ;
           0,      0,    -pc,       ;
          ];
    
H = Hc + Hp + H0;

end

