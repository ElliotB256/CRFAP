function [ H ] = BottomOfShell( t, F, omega0, RFs, gFuBBx, gFuBBy, gFuBBz, py, pz, phase )
%BOTTOMOFSHELL Get H(t) for RF of general polarisation.
%  Evaluates the Hamiltonian for atoms at the specified time. The
%  classical field is of the form:
%   [ Bx sin(wt+phase_i); By sin(wt+py+phase_i); Bz sin(wt+pz+phase_i) ]
% 
%  This function is general for all F. It must be evaluated at the bottom
%  of the shell, ie along z-axis.

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.
swx = sin(RFs.*t+phase);
swy = sin(RFs.*t+py+phase);
swz = sin(RFs.*t+pz+phase);

% Calculate coefficients of Fz, F+, F-
cFz = sum(gFuBBz .* swz, 1);
cFp = sum((gFuBBx .* swx - 1i .* gFuBBy .* swy)/2, 1);
cFm = sum((gFuBBx .* swx + 1i .* gFuBBy .* swy)/2, 1);

% Calculate the Clebsch-Gordon coefficients.
%  Optimisation: Because these coefficients are the same for $\pm$, we only
%  calculate them once.
mFp = -F:F-1;
CGc = ( F .* (F+1) - mFp .* (mFp + 1) ).^0.5;

% Operators for Fz, F+, F-
Fp = diag(CGc, -1);
Fm = diag(CGc, 1);
Fz = diag(-F:1:F);

% Assemble Hamiltonian
H = (omega0 + cFz) * Fz + cFp * Fp + cFm * Fm;

end

