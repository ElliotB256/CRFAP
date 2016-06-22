%% Meshed RFAPs
% Calculate MRF potentials using meshing to save time and resolve finer features
RFs = [3 3.6 4.2 ]';
Rabi = 0.8 * [ 0.5 0.8 1.1 ]';
%RFs = 4.2;
%BRFs = 0.5 ;
Bs=(2.5:0.2:5);
Bs = 3.9:0.05:4.5;
[ F, B ] = MRF.MeshedQuasiEnergies(ZeemanSplit, RFs, Rabi, 'iterations', 5);
plot(B,F,'.');

%%
% Create energy level ladder structure
fundamental = MRF.GetFundamental(RFs);
offsets = fundamental*(-10:1:10)';
lad = repmat(F, length(offsets), 1) + repmat(offsets, 3, size(F, 2)); 
plot(B,lad','.')

%% 
% Locate the minimum of the potential in the stated range
range = [3.9 4.5];
mask = logical(B > range(1) & B < range(2));

% pick our favourite manifold:
fav = 13; 
%plot(B,lad(fav,:),'.') % for testing to check we have right one!

% locate the minimum:
temp = lad(fav, :); temp(~mask) = max(temp(:));
[~,mini] = min(temp); clear temp;
%plot(B(mask), lad(fav, mask), '.'); hold on; plot(B(mini), lad(fav,mini), 'x'); hold off

% Calculate energy distance to other levels from this minimum. This gives
% possible RF transitions truncated to a finite number of photon processes.
spec = abs(lad(:, mini) - lad(fav,mini));

% remove forbidden transitions
% spec = spec(2:3:end);

spec = unique(spec);
plot(B, lad','.');

evapTrans = 16;
spec(evapTrans)