function [ eigF, eigV ] = compute(context, omega0, theta, gamma)
%COMPUTE Computes dressed eigenenergies, see GetDressedEnergies.
%  The arguments are as follows:
%   context: AP.Calculator to use for the evaluation.
%   omega0 : The zeeman splitting due to static field, in MHz
%   theta  : angle, see AP.Calculator
%   phi    : angle, see AP.Calculator

cH = context.GetHamiltonian(omega0, theta, gamma);
periodicity = 2*pi/Floquet.GetFundamental(context.RF);
cU = Floquet.Propagator(cH, periodicity, context.Atom.F);

% Return eigenvectors of propagator if they are requested.
if nargout > 1
    [ eigV, vals ] = eig(cU);
    [ eigF, j] = sort(angle(diag(vals))/periodicity);
    eigV = eigV(:,j);
else
    eigF = sort(angle(eig(cU)))/periodicity;
end

end