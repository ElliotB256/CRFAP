%% Meshed RFAPs
% Calculate MRF potentials using meshing to save time and resolve finer features
RFs = [3 3.6 4.2 ]';
Rabi = 0.8 * [ 0.5 0.5 1.1 ]';
%RFs = 4.2;
%Rabi = 0.5 ;
Bs = (2.5:0.5:5);
Bs = 3.9:0.1:4.5;
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 6);
plot(B,F,'.');

%%
% Create energy level ladder structure
lad = MRF.ladder(RFs, 5, F);

%% 
% Select the different manifolds. Identify the middle most trapped manifold
trapped = lad(3:3:end, :);
flat = lad(2:3:end, :);
antitrapped = lad(1:3:end, :);
fav = ceil(size(trapped, 1)/2);

plot(B, trapped, 'r'); hold on; plot(B, flat, 'Color', [0.5 0.5 0.5]); plot(B, antitrapped, 'Color', [0 0 0]); plot(B, trapped(fav, :), 'r', 'LineWidth', 2); hold off

%% 
% Locate the potential minimum. Accuracy can be increased through iteration
% number for earlier meshing sequence.

[~,mini] = min(trapped(1,:));

%%
% Calculate the energy distance to other levels from this minimum. This gives
% possible RF transitions truncated to a finite number of photon processes.

%spec = (flat - repmat(trapped(fav, :), size(flat, 1), 1));
spec = (flat(:,mini) - trapped(fav, mini));
spec = sort(abs(spec));