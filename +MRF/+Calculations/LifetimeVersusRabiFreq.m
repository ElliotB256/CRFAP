%% Lifetime Versus Rabi Frequency
% Attempts to plot lifetime in the shell versus Rabi frequency, and relate
% this directly to the observed noise spectrum of the RF.

RF42 = 4.2;
BaseRabi42 = 0.357;

RF41 = 4.1;
BaseRabi41 = 0.3696 / 1.3;

RF412 = 4.12;
BaseRabi412 = 0.3562 / 1.3;

RF40 = 4.0;
BaseRabi40 = 0.3455;

BaseRabi = BaseRabi42;
RF = RF42;

Rabis = (0.5:0.05:1.45)*BaseRabi;
ZeemanSplit = (3.5:0.3:5.5)-0.6;

figure(1);
spectra = [];
for Rabi=Rabis
    [spec, debug] = MRF.Spec.Calc(ZeemanSplit, RF, Rabi, 'qdrpGrad', 62, 'ladderN', 8);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '.-', 'Color', [0.5 0.5 0.5]); hold on;
    plot(debug.B(debug.min), debug.trapped(debug.min), '.', 'Color', [0.8 0.2 0.2]); hold on;
end
hold off;

xlabel('Bare Zeeman splitting (MHz)');
ylabel('Dressed state energies (MHz)');
title('Gravitational sag');

figure(2);
plot(Rabis/BaseRabi, spectra', '-', 'Color', [0.8 0.2 0.2]);

%% 
% We now have a list of Rabi frequencies and the spectrum of transitions
% corresponding to them at  approx n\omega \pm 2rabi. Let's focus on just the ones at
% omega \pm Rabi. I'll also assume that the noise is the same either side
% of the dressing rf as it looks to be predominantly a sort of sideband, so
% I'll just pick +Rabi and look at that first.
s = sort(spectra, 1);
s = s(9,:)

% now create interpolating function to map transition freq to Rabi freq
f2rf = @(x) interp1(s,Rabis,x,'linear');
rf2f = @(x) interp1(Rabis, s, x, 'linear');

%%
% Load the measured noise from the scope.

path = 'C:\Data\RFNoise\RFNPK001.csv';
%path = 'C:\Data\RFNoise\RFND000.csv';
fid = fopen(path, 'r');
headerInfo = textscan(fid, '%s %s', 15, 'delimiter', ',');
fclose(fid);

vals = headerInfo{2};

sampleRate = str2double(vals{8});
N = str2double(vals{6});

fdata = csvread(path, 16, 1);
voltage = fdata(:,1);
t = (1:1:N)/sampleRate;

% Perform fourier transform. Focus on noise near to the transitions we picked above.
[f,a] = Util.getFT(t,voltage);
f = f * 1e-6;

i = find(f2rf(f) > min(Rabis(:)) & f2rf(f) < max(Rabis(:)));

plot(f(i), log10(a(i).*conj(a(i))));

%%
% Put it all together. 

% Load from  Y:\Groups\Atomic & Laser\FootGroup\-174 - Enterprise\Shells\multiRF\loading\logbook2016_07_26, first data set.
[~,js] = sort(data(:,1));
plot(data(js,1), data(js,2), '.-', 'Color', [0.3 0.3 0.3]);

% plot rf noise
hold on;
n = log10(a(i).*conj(a(i)));
plot(f2rf(f(i)) / BaseRabi, -n / max(n(:)) * max(data(:,2)), 'r-');

% Direct from pickup away fromm coils
% plot(f2rf(0.8556)*[1 1]/BaseRabi, ylim, 'r-');
% plot(f2rf(0.9480)*[1 1]/BaseRabi, ylim, 'r-');
% plot(f2rf(1.2178)*[1 1]/BaseRabi, ylim, 'g-');
% plot(f2rf(1.2886)*[1 1]/BaseRabi, ylim, 'b-');

plot(f2rf(4.795)*[1 1]/BaseRabi, ylim, 'r-');
plot(f2rf(0.8556)*[1 1]/BaseRabi, ylim, 'r-');
plot(f2rf(4.9998)*[1 1]/BaseRabi, ylim, 'r-');
hold off

legend('Atom Number after 8s', 'pickup noise');
xlabel('Normalised Rabi freq (to 357 kHz)');

%%
% Plot lifetime versus expected resonant frequency...
% Have to manually generate the conversion functions from above code for
% now! :P
[~,i] = sort(data42(:,1));
plot(rf2f42(data42(i,1)*BaseRabi42), data42(i,2), 'x-');

hold on;
[~,i] = sort(data41(:,1));
plot(rf2f41(data41(i,1)*BaseRabi41), data41(i,2), 'x-');
[~,i] = sort(data412(:,1));
plot(rf2f412(data412(i,1)*BaseRabi412), data412(i,2), 'x-');
[~,i] = sort(data40(:,1));
plot(rf2f40(data40(i,1)*BaseRabi40), data40(i,2)*1e6, 'x-');
hold off;

legend('4.2 MHz', '4.1 MHz', '4.12 MHz', '4.0 MHz');

title('Rabi freq. spectroscopy');
ylabel('atom number after 8s');
xlabel('Expected resonance (MHz)');