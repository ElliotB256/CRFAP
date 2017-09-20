% TAAP Trap Frequency Calculator
import Constants.*

% Calculate trap frequencies in the TAAP
BGrad = 1.17 * 200; % Gauss/cm/A * A
RFAmp = 0.6 ./ (0.7 * 2/3); % Gauss
RF = 3.6;
species = 85;

% Based on 5.75 MHz trap bottom during TOP evap of Rb87 at control voltage
% of 2.6V. Converts TOP control voltage to Gauss.
ControlVoltage2TOP = @(voltage) (5.75 ./ 0.7) .* (voltage ./ 2.6);

controlVoltages = linspace(0, 2.6, 100);
results = struct('fx', {}, 'fy', {}, 'fz', {});

for i=1:length(controlVoltages)
    results(i) = getTAAPFrequencies(RF, ControlVoltage2TOP(controlVoltages(i)), RFAmp, BGrad, 'Species', species);
end

f = [[results.fx];[results.fy];[results.fz]]; f(~isreal(f) & f <= 0) = NaN; 
gf = geomean(f, 1);

figure(1); clf
plot(controlVoltages,[results.fx],'-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]); hold on
plot(controlVoltages,[results.fz],'-', 'LineWidth', 2, 'Color', [0.2 0.2 0.8]);
plot(controlVoltages, (abs(real(gf)) > 0).*gf, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2); hold off

xlabel('TOP Control Voltage (V)'); ylabel('Trap Frequency (Hz)');
title({sprintf('Rb%d: TAAP trap frequencies', species), ...
    sprintf('RF=%.2f MHz, Amplitude=%.2f Gauss, Quad=%.2f G/cm^2', ...
    RF, ...
    RFAmp, ...
    BGrad ...
    )});
set(gcf, 'Color', 'w'); box on;


figure(2); clf
plot(ControlVoltage2TOP(controlVoltages),[results.fx],'-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]); hold on
plot(ControlVoltage2TOP(controlVoltages),[results.fz],'-', 'LineWidth', 2, 'Color', [0.2 0.2 0.8]);
plot(ControlVoltage2TOP(controlVoltages), (abs(real(gf)) > 0).*gf, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2); hold off

xlabel('B_{TOP} (Gauss)'); ylabel('Trap Frequency (Hz)');
title({sprintf('Rb%d: TAAP trap frequencies', species), ...
    sprintf('RF=%.2f MHz, Amplitude=%.2f Gauss, Quad=%.2f G/cm^2', ...
    RF, ...
    RFAmp, ...
    BGrad ...
    )});
set(gcf, 'Color', 'w'); box on;