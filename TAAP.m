% Simulate single frequency shell
import Constants.*

zsf = 0.7; %Mhz/Gauss
BGrad = 84; %Gauss/cm
mass = 87; %amu
BRF = 0.5;


% Numerically integrate shell trap for a variety of different TOP fields
N = 50;
BTouch = 2; %G todo programmatically

tms = [];

trapFx = zeros(0,3);
trapFy = zeros(0,3);
trapFz = zeros(0,3);

BTOPs = fliplr(0:BTouch*1.5/N:BTouch*1.5);
for BTOP=BTOPs
    
    % Get top displacement in microns
    xTOP = BTOP/BGrad * 1e4;
    
    % Define potential
    trap = @(a,b,c,t) ShellTrap(...
        a - xTOP .* cos(2*pi*t), b - xTOP .* sin(2*pi*t),c,...
        zsf, BGrad, 1.4, BRF) + ...
        gpe(c, mass);
    
    % Perform time averaging crudely along line x,y=0 to find trap
    % minimum.
    z = 0:-.1:-500;
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
    probeResolution = 20;
    ps = -probeLength:probeLength*2/probeResolution:probeLength;
    zs = zeros(size(ps));
    Nt = 10;
    
    trapFx(end+1,:) = getTrapFreq(ps, timeAverage(@(t) trap(ps,zs,zs+trapMinZ,t),Nt), mass);
    trapFy(end+1,:) = getTrapFreq(ps, timeAverage(@(t) trap(zs,ps,zs+trapMinZ,t),Nt), mass);
    trapFz(end+1,:) = getTrapFreq(ps+trapMinZ, timeAverage(@(t) trap(zs,zs,ps+trapMinZ,t),Nt), mass);
end

plot(BTOPs,trapFx(:,2),'r-'); hold on
plot(BTOPs,trapFy(:,2),'g-');
plot(BTOPs,trapFz(:,2),'b-'); hold off;

xlim([0 TOP4Touch(1.4, BGrad, zsf)]);

%%
 plot(BTOPs,tms); hold on;
 dz = ((resonantEllipsoidWidth(1.4, BGrad, zsf).^2 - (BTOPs/(BGrad*1e-4)).^2).^0.5)/2;
 plot(BTOPs, -dz, ':'); hold off
