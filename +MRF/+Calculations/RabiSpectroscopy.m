%% Plot the lifetime in shell at different Rabi freqs in terms of frequency space.
% Uses data from 26/27 July 2016. Rabi freqs were measured using rf
% spectroscopy.

ZeemanSplit = (3.5:0.3:5.5)-0.6;
RabiRange = (0.3:0.1:1.6);

fprintf('4.2 MHz.\n')
RF42 = 4.2;
BaseRabi42 = 0.357;
spectra42 = [];
for Rabi=RabiRange*BaseRabi42
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF42, Rabi, 'qdrpGrad', 62, 'ladderN', 40);
    spectra42(:,end+1) = spec;
end
[freqs42, atomN42] = MRF.Calculations.rabispectra(data42(:,1)*BaseRabi42, data42(:,2), RabiRange*BaseRabi42, spectra42);
plot(freqs42, atomN42);

fprintf('4.1 MHz.\n')
RF41 = 4.1;
BaseRabi41 = 0.3696 / 1.3;
spectra41 = [];
for Rabi=RabiRange*BaseRabi41
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF41, Rabi, 'qdrpGrad', 62, 'ladderN', 40);
    spectra41(:,end+1) = spec;
end
[freqs41, atomN41] = MRF.Calculations.rabispectra(data41(:,1)*BaseRabi41, data41(:,2), RabiRange*BaseRabi41, spectra41);
hold on; plot(freqs41, atomN41); hold off;

fprintf('4.12 MHz.\n')
RF412 = 4.12;
BaseRabi412 = 0.3562 / 1.3;
spectra412 = [];
for Rabi=RabiRange*BaseRabi412
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF412, Rabi, 'qdrpGrad', 62, 'ladderN', 40);
    spectra412(:,end+1) = spec;
end
[freqs412, atomN412] = MRF.Calculations.rabispectra(data412(:,1)*BaseRabi412, data412(:,2), RabiRange*BaseRabi412, spectra412);
hold on; plot(freqs412, atomN412); hold off;

fprintf('4.00 MHz.\n')
RF40 = 4.0;
BaseRabi40 = 0.3455;
spectra40 = [];
for Rabi=RabiRange*BaseRabi40
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF40, Rabi, 'qdrpGrad', 62, 'ladderN', 40);
    spectra40(:,end+1) = spec;
end
[freqs40, atomN40] = MRF.Calculations.rabispectra(data40(:,1)*BaseRabi40, data40(:,2), RabiRange*BaseRabi40, spectra40);
hold on; plot(freqs40, 1e6*atomN40); hold off;



