%% RF Spec V Rabi Frequency for single RF
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

RF = 3.6;
Rabis = 0.1:0.05:1;
ZeemanSplit = (3.5:0.3:5.5)-0.6;

figure(1);
spectra = [];
for Rabi=Rabis
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF, Rabi, 'qdrpGrad', 62);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '-', 'Color', [0.5 0.5 0.5]); hold on;
    plot(debug.B(debug.min), debug.trapped(debug.min), '.', 'Color', [0.8 0.2 0.2]); hold on;
end
hold off;

xlabel('Bare Zeeman splitting (MHz)');
ylabel('Dressed state energies (MHz)');
title('Gravitational sag');

figure(2);
plot(Rabis, spectra', '-', 'Color', [0.8 0.2 0.2]);

xlabel('Rabi Frequency (MHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq (4 MHz)');

%%
% Use these results to extract the Rabi freq of the dressing RF (including
% Bloch Siegert shifts)

ResonantFreq = 0.46;

% find the resonant freq in the spectra
finer = @(x) interp1(Rabis, spectra(1,:), x, 'pcubic');

trials = 0:0.0001:ResonantFreq;
result = trials(find(finer(trials) > ResonantFreq, 1, 'first'));

fprintf('Rabi frequency of single RF: %.1f kHz\n', result*1e3)

% Measured resonant frequencies for [3, 3.6, 4.2] are [ 0.39 0.46 0.41 ]
% Giving Rabi frequencies of:
% [ 338.7 401.2 357.0]