function [ H ] = F1GeneralBottomOfShell( t, omega0, RFs, gFuBBx, gFuBBy, gFuBBz, py, pz )
%F1CIRCPOL Get H(t) for RF of general polarisation at bottom of the shell.
%  Evaluates the Hamiltonian for F=1 atoms at the specified time.

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

e = exp(1);

H0 = [ 
  -omega0,          0,      0,  ;
        0,          0,      0,  ;
        0,          0, omega0,  ;
     ];

% Calculate coherence terms
cp    = 2^-0.5 * sum( - gFuBBx .* sin(t .* RFs) + 1i *  gFuBBy .* sin(t .* RFs + py), 1);
cn    = 2^-0.5 * sum( - gFuBBx .* sin(t .* RFs) - 1i *  gFuBBy .* sin(t .* RFs + py), 1);

% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above.
Hc = [     ...
         0,     cp,      0,       ;
        cn,      0,     cp,       ;
         0,     cn,      0,       ;
     ];
    
% oscillating component of rf field parallel to local field
% pc = parallel component
pc = sum(gFuBBz .* cos(RFs .* t + pz), 1);
Hp = [ ...
          pc,      0,      0,       ;
           0,      0,      0,       ;
           0,      0,    -pc,       ;
          ];
    
H = Hc + Hp + H0;

end

