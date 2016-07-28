function [ E, B ] = getTrappedLevel(Bs, RFs, Rabi, qdrpGrad)

% Calculate the dressed energy levels and create ladder structure
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 5, 'qdrpGrad', qdrpGrad);
F = MRF.sortEnergies(B,MRF.ladder(RFs, 10, F));

% Select the trapped manifold. To do this, select the manifold which has a
% predominantly positive second derivative.
dB = B(2:end) - B(1:end-1);
d  = (F(:,2:end) - F(:, 1:end-1))./repmat(dB, size(F,1), 1);
dB2 = (dB(1:end-1) + dB(2:end))/2;
d2 = (d(:,2:end) - d(:,1:end-1))./repmat(dB2, size(d,1), 1);

[~,ti] = max(mean(d2, 2));                  % trapped
[~,uti] = min(mean(d2, 2));                 % untrapped
ind = 1:3;
fli = ind(~(ind == ti) & ~(ind == uti));    % flat

E = F(ti, :);

end