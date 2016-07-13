%% Single RF Rabi Freq
% Determines the Rabi frequency for a single RF shell

% 4.2 MHz:
RF = 4.2;
ZeemanSplit = [3.8:0.2:5.0];
Measured = RF+0.390;

% RF = 3;
% ZeemanSplit = [2.8:0.2:3.8];
% Measured = RF+0.392;

% RF = 3.6;
% ZeemanSplit = [3.0:0.3:4.6];
% Measured = RF+0.410;

% ti = 2 for -ve probe, 3 for +ve probe. Check rf spec diagram has right
% transition highlighted in red!
transitionIndex = 3;

Rabis = [0.1:0.05:1];
%Rabis = R+(-1:0.1:1)*0.05;

% Index of the RF transition we are driving. I will plot this in a
% different color at the end just to highlight it.

getTrans = @(spectra) spectra(transitionIndex, :);

figure(1);
spectra = [];
for Rabi=Rabis
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF, Rabi, 'qdrpGrad', 62);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '.-', 'Color', [0.5 0.5 0.5]); hold on;
    plot(debug.B(debug.min), debug.trapped(debug.min), '.', 'Color', [0.8 0.2 0.2]); hold on;
end
hold off;

xlabel('Bare Zeeman splitting (MHz)');
ylabel('Dressed state energies (MHz)');
title('Gravitational sag');

figure(2);
plot(Rabis, spectra', '-', 'Color', [0.5 0.5 0.5]); hold on
plot(xlim, Measured*[1 1], ':b');
plot(Rabis, getTrans(spectra), '--', 'Color', [0.8 0.5 0.5], 'LineWidth', 1.5); hold off

xlabel('Rabi Frequency (MHz)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq (4 MHz)');

%%
% We now find the Rabi freq that agrees with the measured RF spectroscopy
% value. To do this we first create an interpolating function and then
% iterate over this until we get the value we want.

% Now use interp function to find where spec==measured
st = find(getTrans(spectra) - Measured < 0, 1, 'last');
en = find(getTrans(spectra) - Measured > 0, 1, 'first');
dr = (Rabis(en)-Rabis(st))/1000;
vals = Rabis(st):dr:Rabis(en);

f = interp1(Rabis, getTrans(spectra), vals, 'linear');

% find the closest agreement over finer range:
[~,i] = min(abs(f - Measured));
R = vals(i)

figure(2);
hold on; plot([R R], ylim, 'b:'); hold off