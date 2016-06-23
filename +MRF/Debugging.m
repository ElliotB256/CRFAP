rabi100 = [0.391 0.459 0.410]'; % '100%' reference values measured in single-rf traps (MHz)

% Work with series 2 for now
amp30 = 0.5 * rabi100(1);
amp36 = [0.20 0.40 0.60 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.20 1.30 1.40] * rabi100(2); % the final barrier heights
amp42 = 1 * rabi100(3);

qdrpGrad = 62.4511; % Gauss/cm at 20 A

Rabis = [amp30 amp36(2) amp42]';
Bs = 2.5:0.1:4.5;
RF = [3 3.6 4.2]';
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RF, Rabis, 'iterations', 4);
F = MRF.sortEnergies2(B, MRF.ladder(RF, 10, F));

plot(B, F);