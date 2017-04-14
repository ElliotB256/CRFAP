function [ eigF, eigV ] = compute(B, RF, gFuBB, p, periodicity)
%COMPUTE Computes MRF eigenenergies. Used in GetQuasiEnergies.
% B: zeeman splitting, MHz
% RF: dressing frequencies, MHz
% gFuBB: the value of gF * Bohr magneton * field amplitude B. equal to
%        rabi freq for circ polarisations.
% p: input argument struct from GetQuasiEnergies
% periodicity: periodicity of the MRF system.

cH = MRF.Hamiltonian(B, RF, gFuBB, ...
    'F', p.Results.F, ...
    'theta', p.Results.theta, ...
    'phase', p.Results.phase, ...
    'polarisation', p.Results.polarisation);
cU = MRF.Propagator(cH, periodicity, p.Results.F);

% Return eigenvectors of propagator if they are requested.
if nargout > 1
    [ eigV, vals ] = eig(cU);
    [ eigF, j] = sort(angle(diag(vals))/periodicity);
    eigV = eigV(:,j);
else
    eigF = sort(angle(eig(cU)))/periodicity;
end

end