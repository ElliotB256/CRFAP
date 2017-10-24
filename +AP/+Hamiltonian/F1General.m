function [ H ] = F1General( t, omega0, RFs, theta, gamma, gFuBBx, gFuBBy, gFuBBz, py, pz, phase )
%F1GENERAL Get H(t) for RF of general polarisation.
%  Evaluates the Hamiltonian for F=1 atoms at the specified time. The
%  classical field is of the form:
%   [ Bx sin(wt+phase_i); By sin(wt+py+phase_i); Bz sin(wt+pz+phase_i) ]


% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

sTh = sin(theta);
cTh = cos(theta);
sGa = sin(gamma);
cGa = cos(gamma);
swx = sin(RFs.*t+phase);
swy = sin(RFs.*t+py+phase);
swz = sin(RFs.*t+pz+phase);

H0 = [ 
  -omega0,          0,      0,  ;
        0,          0,      0,  ;
        0,          0, omega0,  ;
     ];

% Calculate coherence terms
a = (gFuBBx .* sGa .* sTh .* swx + gFuBBy .* cGa .* swy - gFuBBz .* cTh .* sGa .* swz);
b = (gFuBBx .* cTh .* swx + gFuBBz .* sTh .* swz);
c21 = sum((2^-0.5)*(-1i * a + b),1);
c12 = sum((2^-0.5)*(1i * a + b),1);

Hc = diag([1 1], 1) * c12 + diag([ 1 1 ], -1) * c21;
 
 % Terms parallel to the Hamiltonian
pt = sum(gFuBBx .* cGa .* sTh .* sin(RFs .* t) - gFuBBy .* sGa .* sin(RFs.*t+py) - gFuBBz .* cGa .* cTh .* sin(RFs.*t+pz),1);
Hp = [ pt, 0, 0 ; 0, 0, 0 ; 0, 0, -pt ];

H = Hc + H0 + Hp;

end

