%% MRF Spectroscopy
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

% Important: Check that every minimum is correctly identified!

RFs = [3 3.6 4.2 ]';
Rabi = 0.4 * [ 0.5 0.5 1.1 ]';
BarrierRabi = [0:0.01:1];
ZeemanSplit = [3.7:0.2:4.3];

spectra = [];

for bRabi=BarrierRabi
    Rabi(2) = bRabi;
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF, Rabi, 'ladderN', 20);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '-'); hold on;
end
hold off;

figure(2);
plot(BarrierRabi * 1e3, spectra', '-', 'Color', [0.8 0.2 0.2]);

xlabel('Barrier Rabi Frequency (kHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq');
ylim([0 10]);