function [ H ] = F1CircPolBottomOfShell( t, omega0, RFs, gFuBB, phase )
%F1CIRCPOLBOTTOMOFSHELL Get H(t) at bottom of shell for circ polarised RF.
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

% Calculate coherence terms
cp    = (2.^0.5) * sum( c .* ( e.^(-1i .* (t .* RFs + phase)) ), 1) / 2;
cn    = (2.^0.5) * sum( c .* ( e.^(1i .* (t .* RFs + phase)) ), 1) / 2;

% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above.
Hc = [     ...
         0,     cp,      0,       ;
        cn,      0,     cp,       ;
         0,     cn,      0,       ;
     ];
    
H = Hc + H0;

end

