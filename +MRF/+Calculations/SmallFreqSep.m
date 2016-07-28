%% Small freq sep
% Plot some characterisations of an MRF potential with dF=0.1

RF = [ 4.0 4.1 4.2 ]';
Bs =  3.8:0.025:4.4;
Bs =  3.9:0.025:4.35;
Rabi = [ 134 0.9*144 150 ]' * 1e-3;
%Rabi(2) = 150 * 1e-3;
% Above line will make an extremely flat potential - but note that this
% means the solutions become degenerate and the whole things looks a bit of
% a mess when trying to sort the eigensolutions.

qdrpGrad = 62.5*0.96*3;

[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RF, Rabi, 'iterations', 4, 'qdrpGrad', qdrpGrad);
%F = F+repmat(MRF.gpe(B, qdrpGrad), size(F, 1), 1);
F2 = MRF.sortEnergies(B,MRF.ladder(RF, 10, F));
F2 = F2 + repmat(MRF.gpe(B, qdrpGrad), size(F2, 1), 1);

plot(B, F2, '.-');