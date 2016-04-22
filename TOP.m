%% TOP
% parameters relating to the TOP trap.

gF = Constants.gF; %MHz/Gauss

% Can determine final TOP field magnitude from point we cut to (6 MHz)
finalCut = 6; %MHz
B0 = finalCut / gF

%Find trap frequencies as a function of quadrupole gradient

BGrad = 0:1:400;

%100^2 = cm to m
U = (gF * 1e6 * Constants.hbar * 2 * pi) * BGrad .^ 2 ./ B0 ./ 4 * (100).^2;
trapFreq = (U ./ (87 * Constants.amu)).^0.5 / 2 / pi;
plot(trapFreq);

ztrapFreq = trapFreq * sqrt(8);

plot(BGrad, trapFreq, 'Color', [0.8 0.0 0.0]); hold on;
plot(BGrad, ztrapFreq, 'Color', [0.0 0.0 0.8]); hold off;

title({'Trap frequencies in the TOP Trap',['B0 = ' num2str(B0, '%.1f') ' G (' num2str(finalCut, '%.1f') ' MHz)']});
xlabel('Quadrupole Gradient (G/cm)', 'Interpreter', 'Latex');
ylabel('Trap Frequency (Hz)', 'Interpreter', 'Latex');

legend('radial', 'vertical');

%% Size of condensate
% Use TFA to approximate size of condensate in the TOP trap

geomFreq = (trapFreq .* trapFreq .* ztrapFreq).^(1/3);

% oscillator length in microns
a = (Constants.hbar ./ (87 * Constants.amu * geomFreq * 2 * pi)).^0.5 * 1e6;

N = 3e5; %condensate size of 3 \times 10 ^5 atoms
aScatt = 100 * Constants.bohr * 1e6; %approx for Rb87 in microns
%chemical potential in kHz:
chemPot = 15^(2/5)/2 * (N * aScatt ./ a).^(2/5) .* 2 * pi .* geomFreq / 1e3;

size = (2 * ((chemPot / 1e3) * 2 * pi * Constants.hbar) ./ (Constants.amu * 87 * (2*pi*trapFreq).^2)).^0.5;
size = size * 1e6; %size in microns;
plot(BGrad, size ./ a);

%This size gives reasonable values
sizeB = 1.719 * ( N * aScatt ./ a).^(1/5).*a;
plot(sizeB);