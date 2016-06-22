%% MRF Spectroscopy
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

% Important: Check that every minimum is correctly identified!

RFs = [3 3.6 4.2 ]';
Rabi = [ 0.5 0.5 1.1 ]';
BaseRabi = [ 0.390 0.460 0.410 ]';
Rabi = Rabi .* BaseRabi;
BarrierPct  = 0 : 0.05: 0.4;
BarrierRabi = BarrierPct .* BaseRabi(2);
ZeemanSplit = 3.9:0.1:4.3;
%ZeemanSplit = 3:0.2:4;

spectra = [];

for bRabi=BarrierRabi
    Rabi(2) = bRabi;
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RFs, Rabi, 'ladderN', 40);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '-'); hold on;
end
%hold off;

figure(2);
plot(BarrierRabi * 1e3, spectra', '-', 'Color', [0.8 0.2 0.2]);

xlabel('Barrier Rabi Frequency (kHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq');
ylim([0 10]);