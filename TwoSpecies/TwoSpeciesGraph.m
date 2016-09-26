%% Two Species
% Create a graph showing the potentials we could make for Rb 85 and 87.
% This will trap the 87 in the potential due to 3.5 MHz, and 85 in the
% potential due to 4.0/.1/.2.

n = [ 35 36 37 38 39 40 41 42 ]';
df = 0.1;
RFs = n*df;
Rabi = df * [ 0.7 0 0 0 0 0.7 0.3 0.7 ]';
qdrpGrad = 200;

Bs = 3.4:0.05:4.4;
warning('This may take a while...');
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 6, 'qdrpGrad', qdrpGrad);
%F = MRF.sortEnergies(B, MRF.ladder(RFs, 10, F));
%lad = MRF.ladder(RFs, 3, F);
%glad = lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1);

plot(B, F);


%%
% A very large number of avoided crossings result. Does that look to cause
% a problem for us? Try performing the same calculation using just two
% dressing rfs. For comparison purposes, I've decided to keep the df the
% same - this is in some way a test of how accurate the numerics are.

n = [ 35 36 37 38 39 40 41 42 ]';
df = 0.1;
RFs = n*df;
Rabi = [ 0.2 0 0 0 0 0.0 0.0 0.2 ]';
qdrpGrad = 200;

Bs = 3.4:0.05:4.4;
warning('This may take a while...');
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 6, 'qdrpGrad', qdrpGrad);
F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 10, F));
lad = MRF.ladder(RFs, 3, F2);
glad = lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1);

plot(B, glad);