%% Trap frequency comparison
% For a given rf trap, which has a higher trap freq, F=1 or F=2 of 87/85?

qdrpGrad = 117;
Rabi87 = 0.4*1.1;
RF     = 4.2;
%Bs=(2.6:0.1:3.4);
Bs=(4.0:0.1:4.5);

[ F1, B1 ] = MRF.MeshedQuasiEnergies(Bs, RF, Rabi87, 'iterations', 6, 'F', 1);
[ F2, B2 ] = MRF.MeshedQuasiEnergies(Bs, RF, 2/3*Rabi87, 'iterations', 6, 'F', 2);

z87 = B1 / qdrpGrad / 0.7 * 1e4;
z85 = B2 / qdrpGrad / (0.7*2/3) * 1e4;

plot(z87-mean(z87), F1, 'r'); hold on;
plot(z85-mean(z85), F2, 'b'); hold off;

trapFreq = Util.getTrapFreq(z87, F1(3,:), 87);
fprintf('Rb87 Trap freq = %.1f.\n', trapFreq(2))

trapFreq1 = Util.getTrapFreq(z85, F2(4,:), 85);
trapFreq2 = Util.getTrapFreq(z85, F2(5,:), 85);
fprintf('Rb85 Trap freq = (%.1f, %.1f).\n', trapFreq1(2), trapFreq2(2))