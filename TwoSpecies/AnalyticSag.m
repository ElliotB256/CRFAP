%% Analytic Sag
% Attempt to determine the conditions to contact two species. Calculations
% from analytic formulae.

deltaF = 0.6;
gF87 = 0.7;
gF85 = 2/3 * gF87;
g = 9.81;

% All distances should be um. Qdrp Grad in G/cm. Magnetic fields in Gauss.

% Get shifts in the resonance position due to other field.
resPosShift87 = @(B1, qdrpGrad) gF87 * B1 .^ 2 / ( 2 * deltaF * qdrpGrad) * 1e4;
resPosShift85 = @(B2, qdrpGrad) gF85 * B2 .^ 2 / ( 2 * deltaF * qdrpGrad) * 1e4;

% Get trap frequencies in single rf ap
trap85 = @(B1, qdrpGrad) qdrpGrad .^ 2 * gF85 * Constants.h / (85 * Constants.amu) / B1 * 1e4;
sag85 = @(B1, qdrpGrad) g / ( 2 * qdrpGrad .^ 2 * gF85 * 2 / (Constants.amu * 85 * B1) );