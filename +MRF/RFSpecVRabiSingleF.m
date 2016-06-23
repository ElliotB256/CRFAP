%% RF Spec V Rabi Frequency for single RF
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

RF = 4;
Rabis = [0.1:0.1:3];
ZeemanSplit = [3.5:0.3:5.5];

figure(1);
spectra = [];
for Rabi=Rabis
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF, Rabi, 'qdrpGrad', 62);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '-', 'Color', [0.5 0.5 0.5]); hold on;
    plot(debug.B(debug.min), debug.trapped(debug.min), '.', 'Color', [0.8 0.2 0.2]); hold on;
end
hold off;

figure(2);
plot(Rabis, spectra', '-', 'Color', [0.8 0.2 0.2]);

xlabel('Rabi Frequency (MHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq (4 MHz)');

%%
% Use these results to extract the Rabi freq of the dressing RF (including
% Bloch Siegert shifts)