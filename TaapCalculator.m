% TAAP Trap Frequency Calculator
import Constants.*

% Species: 85 or 87
species = 87;
%RB85
BGrad = 84;%;1.17 * 100; % Gauss/cm/A * A
RFAmp = 0.5;%0.250 ./ 0.7 * 3/2; % Gauss
RF = 1.4%3.0;

switch species
    case 85
        zsf = 0.7 * 2/3;
        mFtilde = 2;
        mass = 85; %amu
    case 87
        zsf = 0.7;
        mFtilde = 1;
        mass = 87; %amu
    otherwise
        error('Unknown species');
end

% Based on 5.75 MHz trap bottom during TOP evap. Converts control voltage
% to TOP field in Gauss
ControlVoltage2TOP = @(voltage) 5.75 ./ 0.7 .* (voltage ./ 2.6);

controlVoltages = linspace(0, 2.6, 100);
trapFx = zeros(length(controlVoltages),3);
trapFy = zeros(length(controlVoltages),3);
trapFz = zeros(length(controlVoltages),3);
tms = zeros(length(controlVoltages),0);

for i=1:length(controlVoltages)
    BTOP = ControlVoltage2TOP(controlVoltages(i));
    
    % Get displacement of quadrupole centre in microns
    xTOP = BTOP/BGrad * 1e4;
    
    % Define potential
    trap = @(a,b,c,t) SRF.ShellTrap(...
        a - xTOP .* cos(2*pi*t), b - xTOP .* sin(2*pi*t),c,...
        zsf, BGrad, RF, RFAmp, mFtilde) + ...
        gpe(c, mass);
    
    % Perform time averaging crudely along line x,y=0 to find trap
    % minimum.
    z = 0:-.1:-2000;
    x = zeros(size(z));
    y = x;
    
    pot = timeAverage(@(t) trap(x,y,z,t), 10);
    
    [~,ip] = min(pot);
    trapMinZ = z(ip);
    tms(i) = trapMinZ;
    
    clear ip pot;
    
    % Time average along x,y,z axis at the minimum spot.
    % Take a small size around the trap center.
    probeLength = 3; %um
    probeResolution = 100;
    ps = -probeLength:probeLength*2/probeResolution:probeLength;
    zs = zeros(size(ps));
    Nt = 10;
    
    trapFx(i,:) = Util.getTrapFreq(ps, timeAverage(@(t) trap(ps,zs,zs+trapMinZ,t),Nt), mass);
    trapFy(i,:) = Util.getTrapFreq(ps, timeAverage(@(t) trap(zs,ps,zs+trapMinZ,t),Nt), mass);
    trapFz(i,:) = Util.getTrapFreq(ps+trapMinZ, timeAverage(@(t) trap(zs,zs,ps+trapMinZ,t),Nt), mass);
end

figure(1); clf
plot(controlVoltages,trapFx(:,2),'-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]); hold on
plot(controlVoltages,trapFz(:,2),'-', 'LineWidth', 2, 'Color', [0.2 0.2 0.8]);

%geometric frequency
plot(controlVoltages, (abs(real(trapFx(:,2))) > 0).*(trapFx(:,2) .* trapFy(:,2) .* trapFz(:,2)).^(1/3), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
 hold off
xlabel('TOP Control Voltage (V)');
ylabel('Trap Frequency (Hz)');
title({sprintf('Rb%d: TAAP trap frequencies', species), ...
    sprintf('RF %.2f MHz, Rabi freq %.0f kHz, Quad %.2f G/cm^2', ...
    RF, ...
    RFAmp * zsf * 1000, ...
    BGrad ...
    )});

figure(2); clf;
plot(ControlVoltage2TOP(controlVoltages),trapFx(:,2),'-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]); hold on
plot(ControlVoltage2TOP(controlVoltages),trapFz(:,2),'-', 'LineWidth', 2, 'Color', [0.2 0.2 0.8]);

%geometric frequency
plot(ControlVoltage2TOP(controlVoltages), (abs(real(trapFx(:,2))) > 0).*(trapFx(:,2) .* trapFy(:,2) .* trapFz(:,2)).^(1/3), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
 hold off
xlabel('TOP Field (Gauss)');
ylabel('Trap Frequency (Hz)');
title({sprintf('Rb%d: TAAP trap frequencies', species), ...
    sprintf('RF %.2f MHz, Rabi freq %.0f kHz, Quad %.2f G/cm^2', ...
    RF, ...
    RFAmp * zsf * 1000, ...
    BGrad ...
    )});