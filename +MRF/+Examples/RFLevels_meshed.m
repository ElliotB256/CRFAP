%% Meshed RFAPs
% Calculate MRF potentials using meshing to save time and resolve finer features
RFs = [3 3.6 4.2 ]';
Rabi = 0.3 * [ 1 0.7 1 ]';
qdrpGrad = 100;

Bs=(2.5:0.2:5);
%Bs=(0:0.2:7);
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 3, 'qdrpGrad', qdrpGrad, 'F', 1);

%%
% Make a ladder, sort energies, then re-ladder to get correct states. Add
% in gpe.
F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 30, F));
lad = MRF.ladder(RFs, 30, F2);
%plot(B, lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1), '-');
plot(B, lad, '-', 'Color', [0.7 0.7 0.7]);
xlim([2 5]);
ylim([14.2 15.8]);
hold on; plot(B, lad(33,:), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2); hold off;
set(gca, 'YTick', 14.4 + [0 0.6 1.2 ]);
set(gca, 'YTickLabel', { '0', '0.6', '1.2' } );
xlabel('Zeeman splitting (MHz)', 'FontSize', 14);
ylabel('V/h (MHz)', 'FontSize', 14);