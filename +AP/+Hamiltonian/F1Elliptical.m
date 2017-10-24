function [ H ] = F1Elliptical( t, omega0, RFs, theta, gamma, gFuBBx, gFuBBy )
%F1ELLIPTICAL Get H(t) for RF of elliptical polarisation.
%  Evaluates the Hamiltonian for F=1 atoms at the specified time.
%  Assumes:
%   -Bz=0
%   -phi_y = pi/2, such that the dressing RF is more like an 'unbalanced'
%    circular polarised case.

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

cTh = cos(theta);
sTh = sin(theta);
cGa = cos(gamma);
sGa = sin(gamma);

H0 = [ 
  -omega0,          0,      0,  ;
        0,          0,      0,  ;
        0,          0, omega0,  ;
     ];

% Calculate coherence terms
cp = 2^-0.5 * sum(1i * gFuBBy .* cGa .* cos(t .* RFs) + gFuBBx .* ( cTh + 1i .* sGa .* sTh ) .* cos(t .* RFs),1);
cn = 2^-0.5 * sum(-1i * gFuBBy .* cGa .* cos(t .* RFs) + gFuBBx .* ( cTh - 1i .* sGa .* sTh ) .* cos(t .* RFs),1);

% Note: The Clebsch-Gordon coefficient, here sqrt(2), is incorporated
% above.
Hc = [     ...
         0,     cp,      0,       ;
        cn,      0,     cp,       ;
         0,     cn,      0,       ;
     ];
    
% oscillating component of rf field parallel to local field
% pc = parallel component
pc = sum( - gFuBBy .* sGa .* cos(RFs .* t) + gFuBBx .* sin(RFs .* t) .* cGa .* sTh, 1);

Hp = [ ...
          pc,      0,      0,       ;
           0,      0,      0,       ;
           0,      0,    -pc,       ;
          ];
    
H = Hc + Hp + H0;

end

