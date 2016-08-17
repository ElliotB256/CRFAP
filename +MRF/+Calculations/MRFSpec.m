%% MRF Spectroscopy
% Calculates the RF spectroscopy versus rabi frequency for a single RF. See
% this as more of a test that the functions are working correctly.

%%
% Load the most recent quadrupole gradient calibration and rabiFreqs.

thisLoc = mfilename('fullpath');

quadCalFile = [fileparts(thisLoc) filesep 'quadGradient.mat'];
if ~exist(quadCalFile, 'file')
    error('run QuadGradient_top.m first!');
end
load(quadCalFile);

rabiFreqFile = [fileparts(thisLoc) filesep 'rabiFreqs.mat'];
if ~exist(rabiFreqFile, 'file')
    error('run SingleRFRabiFreqs.m first!');
end
load(rabiFreqFile);

clear rabiFreqFile quadCalFile

Rabi = Rabi';
%MRFSpecQuad = MRFSpecQuad * 0.65;

%%

RFs = [ 3 3.6 4.2 ]';
AmplitudeMod = [ 0.5*1.2 0.5 1.0 ]';

% TODO: more exact implementation that scales with barrier height
amplifierFactor = 0.98*[1 1 1]'; %[ 0.98 0.98 0.98 ]';

ActualRabi = AmplitudeMod .* Rabi .* amplifierFactor;
BarrierPct  = 0 : 0.05: 0.65;
BarrierRabi = BarrierPct .* Rabi(2);
ZeemanSplit = 3.8:0.1:4.5;

spectra = [];

for bRabi=BarrierRabi
    fprintf('Calculating bRabi=%.0f kHz\n', bRabi * 1000)
    ActualRabi(2) = bRabi * amplifierFactor(2);
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RFs, ActualRabi, 'ladderN', 40, 'qdrpGrad', MRFSpecQuad, 'iterations', 10);
    spectra(:,end+1) = spec;

    plot(debug.B, debug.trapped, '.-', 'Color', [0.5 0.5 0.5]); hold on;
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