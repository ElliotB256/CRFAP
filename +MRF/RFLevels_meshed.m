%% Meshed RFAPs
% Calculate MRF potentials using meshing to save time and resolve finer features
RFs = [3 3.6 4.2 ]';
Rabi = 0.8 * [ 0.5 0 1.1 ]';
qdrpGrad = 100;

Bs=(2.5:0.2:5);
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 4, 'qdrpGrad', qdrpGrad);

%%
% Make a ladder, sort energies, then re-ladder to get correct states. Add
% in gpe.
F = MRF.sortEnergies2(B, MRF.ladder(RFs, 3, F));
lad = MRF.ladder(RFs, 3, F);
plot(B, lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1), '-');