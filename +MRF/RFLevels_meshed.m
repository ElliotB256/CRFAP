%%
% With more general propagator notation
RFs = [3 3.6 4.2 ]';
BRFs = 0.8 * [ 0.5 0.2 1.1 ]';
periodicity = 10*pi/3;
deltaB = 0.1;
Bs=(2.5:deltaB:5);

eigF2 = MRF.GetQuasiEnergies(Bs, RFs, BRFs);

iterations = 2;
deltaB = Bs(2)-Bs(1);

for i=1:iterations
deltaB = deltaB / 2;
% locate interesting points to mesh around. Here we find stationary points.
d = diff(eigF2, 1, 2)';
d2 = diff(sign(d), 1, 1);
j = find(abs(d2(:,1)) > 0.5);
j2 = min(j+1, length(Bs));

% mesh around the stationary points
mesh = [deltaB deltaB];
cB = unique([Bs(j) Bs(j2)]);
Bs2 = repmat(cB, size(mesh, 2), 1) + repmat(mesh', 1, size(cB, 2));
Bs2 = Bs2(:)';

% Calculate energies of these new points and add to the list.
eigF2 = [eigF2 MRF.GetQuasiEnergies(Bs2, RFs, BRFs)];
Bs = [Bs Bs2];

end

plot(Bs, eigF2', '.')

%%
Bs=(2.5:0.2:5);
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, BRFs, 'iterations', 4);
plot(B,F,'.');