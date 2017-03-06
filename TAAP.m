% Simulate single frequency shell
import Constants.*

zsf = Constants.zeemansplit; %Mhz/Gauss
BGrad = 1.1731 * 100; %2.42*100; %Gauss/cm
mass = 87; %amu
RFAmp = 0.400 ./ zsf;
RF = 4.2;

VITopVoltage = 0.6;
VI2TOP = @(vi) 5.75 ./ zsf ./ 2.6 .* vi;

% Numerically integrate shell trap for a variety of different TOP fields
N = 150;
BTouch = SRF.TOP4Touch(RF, BGrad, zsf); %G todo programmatically

tms = [];

trapFx = zeros(0,3);
trapFy = zeros(0,3);
trapFz = zeros(0,3);

%Bend = 1;%BTouch*1.5;
Bend = 1.5*SRF.TOP4Touch(RF, BGrad, zsf);

BTOPs = fliplr(0:Bend/N:Bend);
for BTOP=BTOPs
    
    % Get top displacement in microns
    xTOP = BTOP/BGrad * 1e4;
    
    % Define potential
    trap = @(a,b,c,t) SRF.ShellTrap(...
        a - xTOP .* cos(2*pi*t), b - xTOP .* sin(2*pi*t),c,...
        zsf, BGrad, RF, RFAmp) + ...
        gpe(c, mass);
    
    % Perform time averaging crudely along line x,y=0 to find trap
    % minimum.
    z = 0:-.1:-2000;
    x = zeros(size(z));
    y = x;
    
    pot = timeAverage(@(t) trap(x,y,z,t), 10);
    
    [~,ip] = min(pot);
    trapMinZ = z(ip);
    tms(end+1) = trapMinZ;
    
    clear ip pot;
    
    % Time average along x,y,z axis at the minimum spot.
    % Take a small size around the trap center.
    probeLength = 3; %um
    probeResolution = 100;
    ps = -probeLength:probeLength*2/probeResolution:probeLength;
    zs = zeros(size(ps));
    Nt = 10;
    
    trapFx(end+1,:) = Util.getTrapFreq(ps, timeAverage(@(t) trap(ps,zs,zs+trapMinZ,t),Nt), mass);
    trapFy(end+1,:) = Util.getTrapFreq(ps, timeAverage(@(t) trap(zs,ps,zs+trapMinZ,t),Nt), mass);
    trapFz(end+1,:) = Util.getTrapFreq(ps+trapMinZ, timeAverage(@(t) trap(zs,zs,ps+trapMinZ,t),Nt), mass);
end

plot(BTOPs,trapFx(:,2),'-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]); hold on
plot(BTOPs,trapFz(:,2),'-', 'LineWidth', 2, 'Color', [0.2 0.2 0.8]);
plot([1 1] * VI2TOP(0.6), ylim, 'k-');

%geometric frequency
plot(BTOPs, (abs(real(trapFx(:,2))) > 0).*(trapFx(:,2) .* trapFy(:,2) .* trapFz(:,2)).^(1/3), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
 hold off


xlim([0 Bend]);

title({'TAAP trap frequencies', ...
    sprintf('RF %.2f MHz, Rabi freq %.0f kHz, Quad %.2f G/cm^2', ...
    RF, ...
    RFAmp * Constants.zeemansplit * 1000, ...
    BGrad ...
    )});

%%
%  plot(BTOPs,tms); hold on;
%  dz = ((resonantEllipsoidWidth(2, BGrad, zsf).^2 - (BTOPs/(BGrad*1e-4)).^2).^0.5)/2;
%  plot(BTOPs, -dz, ':'); hold off
