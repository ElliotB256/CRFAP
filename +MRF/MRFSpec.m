%% MRF Spectroscopy
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

% Important: Check that every minimum is correctly identified!

RFs = [3 3.6 4.2 ]';
Rabi = [ 0.5 0.5 1.1 ]';
BaseRabi = [ 0.390 0.460 0.410 ]';
BaseRabi = [ 0.341 0.357 0.357 ]';
Rabi = Rabi .* BaseRabi;
BarrierPct  = 0 : 0.05: 0.65;
BarrierRabi = BarrierPct .* BaseRabi(2);
ZeemanSplit = 3.8:0.1:4.5;
%ZeemanSplit = 3:0.2:4;
QdrpGrad = 62.5;

spectra = [];

for bRabi=BarrierRabi
    Rabi(2) = bRabi;
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RFs, Rabi, 'ladderN', 40, 'qdrpGrad', QdrpGrad);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '-', 'Color', [0.5 0.5 0.5]); hold on;
    plot(debug.B(debug.min), debug.trapped(debug.min), '.', 'Color', [0.8 0.2 0.2]); hold on;
end
hold off;
title('MRF APs used for rf spectroscopy calc');
xlabel('Bare Zeeman Splitting (MHz)');
ylabel('Energy (MHz)');

figure(2);
plot(BarrierRabi * 1e3, spectra', '.-', 'Color', [0.8 0.2 0.2]); ylim([3.5 5.5]);
ylim([4.5 4.7]);

xlabel('Barrier Rabi Frequency (kHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq');