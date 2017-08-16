function [ H ] = F1LinPol( t, omega0, RFs, theta, gamma, gFuBB, phase )
%F1LINPOL Get H(t) for lin polarised RF, F=1.
% 
% gFuBB: the amplitude of the rf field expressed in MHz

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

% precalculate trigonometric quantities
c = gFuBB;
sTh = sin(theta);
cTh = cos(theta);
sGa = sin(gamma);
cGa = cos(gamma);

H0 = [ 
       omega0,      0,      0,       ;
            0,      0,      0,       ;
            0,      0,-omega0,       ;
     ];

% Coherence terms:
% (cos phi + i sin phi) * cos theta * gFuBB / sqrt(2) * sin(wt+phase)
pc = sum( (cTh + 1i .* sGa .* sTh) .* sin(t .* RFs + phase)) .* c / (2.^0.5);
nc = sum( (cTh - 1i .* sGa .* sTh) .* sin(t .* RFs + phase)) .* c / (2.^0.5);

Hc = [     ...
         0,     pc,      0,       ;
        nc,      0,     pc,       ;
         0,     nc,      0,       ;
     ];
 
% Component of field parallel to quantisation axis:
hz = sum(c .* cGa .* sTh .* sin(RFs .* t + phase), 1);
Hp = [ ...
          hz,      0,      0,       ;
           0,      0,      0,       ;
           0,      0,    -hz,       ;
          ];
    
H = Hc + H0 + Hp;

end