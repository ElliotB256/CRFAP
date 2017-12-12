function Sample(instance)
%SAMPLE Calculate eigenenergies at each B.

theta = instance.Theta;
gamma = instance.Gamma;
B = repmat(instance.RF, length(theta), 1);

if (instance.Verbose)
    fprintf('Calculating energies at %d field points.\n', length(theta))
end

% Calculate first set of eigenvalues at initial points
[eigE,eigV] = instance.APCalculator.GetDressedEnergies(B, theta, gamma);

instance.mE = eigE;
instance.mEV = eigV;
instance.Dirty = 0;

end