function Sample(instance)
%SAMPLE Calculate eigenenergies at each B.

B = instance.StartB;

if (instance.Verbose)
    fprintf('Calculating initial energies at %d field points.\n', length(B))
end
% Calculate first set of eigenvalues at initial points
[eigE,eigV] = instance.APCalculator.GetDressedEnergies(B, 0, 0);

% Assume: B is uniformly spaced.
deltaB = instance.StartB(2)-instance.StartB(1);

for i=1:instance.MeshIterations
deltaB = deltaB / 2;

cB = Util.Refine(B, eigE(1,:));

if instance.QuadGrad > 0.0001
    % if Quadrupole gradient is specified then also mesh according to gravity.
    % convert Zeeman splitting into microns, then calculate gpe.
    gp = MRF.gpe(B, instance.QuadGrad, 'gF', abs(instance.APCalculator.Atom.gFuB), 'mass', instance.APCalculator.Atom.Mass);
    
    % Modify eigen energies to account for energy shift
    modEig = eigE + repmat(gp, instance.APCalculator.HilbertSpaceSize, 1);
    
    % Collect unique values suggested from eigenvalues
    gcB = [];
    for j=1:size(eigE, 1)
        gcB = [gcB Util.Refine(B, modEig(j,:))];
    end
    gcB = uniquetol(gcB, deltaB/20);
    
    % Remove duplicate elements
    cB = Util.UniquePick( B, [ cB gcB ]);
end

% Calculate energies of these new points and add to the list.
Bs2 = cB(:)';
if (instance.Verbose)
    fprintf('Calculating meshed energies at %d field points.\n', length(Bs2))
end
[eigF2,eigV2] = instance.APCalculator.GetDressedEnergies(Bs2, 0, 0);
eigE = [eigE eigF2];
B = [B Bs2];

% sort results by B
[~,j] = sort(B);
B = B(j);
eigE = eigE(:,j);

end

%Store results in sampler object.
instance.mB = B;
instance.mE = eigE;
instance.Dirty = 0;

end