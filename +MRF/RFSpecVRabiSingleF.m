%% RF Spec V Rabi Frequency for single RF
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

RF = 4.2;
Rabis = [0.1:0.1:1];
ZeemanSplit = [3.5:0.2:4.5];

spectra = [];
for Rabi=Rabis
spectra(:,end+1) = MRF.Spec.Calc(ZeemanSplit, RF, Rabi);
end

plot(Rabis, spectra', '-', 'Color', [0.8 0.2 0.2]);

xlabel('Rabi Frequency (MHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq');