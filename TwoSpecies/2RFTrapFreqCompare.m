%% Trap frequency comparison
% For a given rf trap, which has a higher trap freq, F=1 or F=2 of 87/85?

qdrpGrad = 100.0;
Rabi87  = [ 0.14 0.4 ]'; 
RFs     = [ 3   4.5  ]';
Bs87=(4.2:0.05:4.7);
Bs85=(2.9:0.05:3.2);

[ F1, B1 ] = MRF.MeshedQuasiEnergies(Bs87, RFs, Rabi87, 'iterations', 6, 'F', 1, 'QdrpGrad', qdrpGrad, 'gF', 0.7);
[ F2, B2 ] = MRF.MeshedQuasiEnergies(Bs85, RFs, 2/3*Rabi87, 'iterations', 6, 'F', 2, 'QdrpGrad', qdrpGrad, 'gF', 0.7*2/3);

z87 = B1 / qdrpGrad / 0.7 * 1e4;
z85 = B2 / qdrpGrad / (0.7*2/3) * 1e4;

plot(z87, F1, 'r'); hold on;
plot(z85, F2, 'b'); hold off;

trapFreq = Util.getTrapFreq(z87, F1(3,:), 87);
fprintf('Rb87 Trap freq = %.1f.\n', trapFreq(2))

trapFreq1 = Util.getTrapFreq(z85, F2(4,:), 85);
trapFreq2 = Util.getTrapFreq(z85, F2(5,:), 85);
fprintf('Rb85 Trap freq = (%.1f, %.1f).\n', trapFreq1(2), trapFreq2(2))

%
% Plot graphs including gravitational shift
F1g  = F1 + repmat(MRF.gpe(B1, qdrpGrad, 'gF', 0.7), 3, 1);
F2g  = F2 + repmat(MRF.gpe(B2, qdrpGrad, 'gF', 0.7 * 2 / 3), 5, 1);

plot(z87, F1g, 'r.-'); hold on;
plot(z85, F2g, 'b.-'); hold on;
[~,i] = min(F1g(3,:)); plot(z87(i), F1g(3,i), 'ro', 'MarkerSize', 12); hold on;
[~,i] = min(F2g(5,:)); plot(z85(i), F2g(5,i), 'bo', 'MarkerSize', 12); hold on;