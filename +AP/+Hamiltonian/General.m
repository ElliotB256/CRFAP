function [ H ] = General( t, F, omega0, RFs, theta, gamma, gFuBBx, gFuBBy, gFuBBz, py, pz, phase )
%GENERAL Get H(t) for RF of general polarisation.
%  Evaluates the Hamiltonian for atoms at the specified time. The
%  classical field is of the form:
%   [ Bx sin(wt+phase_i); By sin(wt+py+phase_i); Bz sin(wt+pz+phase_i) ]
%  
%  The angles theta and gamma are defined in the sherical quadrupole such
%  that the coordinates (x,y,z) are:
%     x = d sin(theta) cos(gamma)
%     y = d sin(theta) sin(gamma)
%     z = d cos(theta) / 2
%    For 'spheroid radius' d. The quadrupole field is
%     Bquad = B' (x, y, -2z)
%    In this picture the quantisation axis at a point p=(x,y,z) is:
%      [ cos(gamma) sin(theta), sin(gamma) sin(theta), -cos(theta) ]
%  
%   At the bottom of the shell we have theta=pi and gamma=0.
%  
%  This function is general for all F.

% THIS FUNCTION WILL BE CALLED FREQUENTLY. Make it as lean as possible.

sTh = sin(theta);
cTh = cos(theta);
sGa = sin(gamma);
cGa = cos(gamma);
swx = sin(RFs.*t+phase);
swy = sin(RFs.*t+py+phase);
swz = sin(RFs.*t+pz+phase);

Tx = -swx .* gFuBBx;
Ty = -swy .* gFuBBy;
Tz = -swz .* gFuBBz;

% Calculate coefficients of Fz, F+, F-
cFz = sum(-gFuBBx .* cGa .* sTh .* swx + gFuBBy .* sGa .* swy + gFuBBz .* cGa .* cTh .* swz, 1);
cFp = sum(gFuBBx .* cTh .* swx / 2 - 1i / 2 .* gFuBBx .* sGa .* sTh .* swx ... 
            - 1i / 2 .* gFuBBy .* cGa .* swy + ...
            1i / 2 * gFuBBz .* cTh .* sGa .* swz + 1/2 * gFuBBz .* sTh .* swz, 1);
cFm = sum(gFuBBx .* cTh .* swx / 2 + 1i / 2 .* gFuBBx .* sGa .* sTh .* swx ... 
            + 1i / 2 .* gFuBBy .* cGa .* swy ...
            - 1i / 2 * gFuBBz .* cTh .* sGa .* swz + 1/2 * gFuBBz .* sTh .* swz, 1);

cFz = sum( ...
    Tx .* cGa .* sTh + ...
    Ty .* sGa .* sTh + ...
    Tz        .* cTh  ...
    , 1);

cFp = sum( ...
    Tx .* (cGa .* cTh - 1i .* sGa)/2 + ...
    Ty .* (sGa .* cTh + 1i .* cGa)/2 + ...
    Tz .* (sTh)/2 ...
    , 1);
        
cFm = sum( ...
    Tx .* (cGa .* cTh + 1i .* sGa)/2 + ...
    Ty .* (sGa .* cTh - 1i .* cGa)/2 + ...
    Tz .* (+sTh)/2 ...
    , 1);

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

