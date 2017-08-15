function [ H ] = F1CircPol( t, omega0, RFs, gFuBB, theta, phase )
%F1CIRCPOL Get H(t) for circ polarised RF.
%  Evaluates the Hamiltonian for F=1 atoms at the specified time. All MRF
%  components are circularly polarised.

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

e = exp(1);
c = gFuBB;

H0 = [ 
   omega0,          0,      0,  ;
        0,          0,      0,  ;
        0,          0,-omega0,  ;
     ];
         
% Assume: All rf components are circ pol. \alpha_i, \beta_i \to \alpha, \beta
alpha = (cos(theta) - 1)/2;
beta  = (cos(theta) + 1)/2;

% Calculate coherence terms
cp    = (2.^0.5) * sum( c .* ( alpha .* e.^(1i .* (t .* RFs + phase)) + beta .* e.^(-1i .* (t .* RFs + phase)) ), 1) / 2;
cn    = (2.^0.5) * sum( c .* ( alpha .* e.^(-1i .* (t .* RFs + phase)) + beta .* e.^(1i .* (t .* RFs + phase)) ), 1) / 2;

% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above.
Hc = [     ...
         0,     cp,      0,       ;
        cn,      0,     cp,       ;
         0,     cn,      0,       ;
     ];
    
% oscillating component of rf field parallel to local field
% pc = parallel component
pc = sum(c .* sin(theta) .* cos(RFs .* t + phase), 1);
Hp = [ ...
          pc,      0,      0,       ;
           0,      0,      0,       ;
           0,      0,    -pc,       ;
          ];
    
H = Hc + Hp + H0;

end

