function [ H ] = F1EllipticalBottomOfShell( t, omega0, RFs, gFuBBx, gFuBBy )
%F1ELLIPTICALBOTTOMOFSHELL Get H(t) for RF of elliptical polarisation.
%  Evaluates the Hamiltonian for F=1 atoms at the specified time.
%  Assumes:
%   -Bz=0
%   -phi_y = pi/2, such that the dressing RF is more like an 'unbalanced'
%    circular polarised case.
%   -theta,phi -> 0, ie bottom of shell.

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

H0 = [ 
  -omega0,          0,      0,  ;
        0,          0,      0,  ;
        0,          0, omega0,  ;
     ];

% Calculate coherence terms
cp = 2^-0.5 * sum(1i * gFuBBy .* cos(t .* RFs) - gFuBBx .* sin(t .* RFs),1);
cn = 2^-0.5 * sum(-1i * gFuBBy .* cos(t .* RFs) - gFuBBx .* sin(t .* RFs),1);

% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above.
Hc = [     ...
         0,     cp,      0,       ;
        cn,      0,     cp,       ;
         0,     cn,      0,       ;
     ];
    
H = Hc + H0;

end

