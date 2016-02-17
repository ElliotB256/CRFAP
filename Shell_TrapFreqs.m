% Simulate single frequency shell
import Constants.*

zsf = Constants.zeemansplit; %Mhz/Gauss
BGrad = 84; %Gauss/cm
rb = rubidium87();
mass = rb.mass; %amu
BRF = 0.5;
RF = 4;


tms = [];

trapFx = zeros(0,3);
trapFy = zeros(0,3);
trapFz = zeros(0,3);

BTOP = 0;
   
BGrads = 30:5:300;

for BGrad=BGrads
    
    % Get top displacement in microns
    xTOP = BTOP/BGrad * 1e4;
    
    % Define potential
    trap = @(a,b,c) ShellTrap(...
        a, b, c,...
        zsf, BGrad, RF, BRF) + ...
        gpe(c, mass);
    
    % first iter: crudely find trap minimum.
    z = 0:-1:-1000;
    x = zeros(size(z)); y = x;
    pot = trap(x,y,z);
    
    [~,ip] = min(pot);
    trapMinZ = z(ip);

    % second iter: refine trap minimum
    z = (-3:0.001:3) + trapMinZ;
    x = zeros(size(z)); y = x;
    pot = trap(x,y,z);
    [~,ip] = min(pot);
    trapMinZ = z(ip);
    tms(end+1) = trapMinZ;
    
    clear ip pot;
    
    % Time average along x,y,z axis at the minimum spot.
    % Take a small size around the trap center.
    probeLength = 5; %um
    probeResolution = 100;
    ps = -probeLength:probeLength*2/probeResolution:probeLength;
    zs = zeros(size(ps));
    Nt = 10;
    
    trapFx(end+1,:) = getTrapFreq(ps,trap(ps,zs,zs+trapMinZ), mass);
    trapFy(end+1,:) = getTrapFreq(ps,trap(ps,zs,zs+trapMinZ), mass);
    trapFz(end+1,:) = getTrapFreq(ps+trapMinZ,trap(zs,zs,ps+trapMinZ), mass);
    
end

%%

figure;
actualH = plot(BGrads, tms); hold on
theoryH = plot(BGrads, -1e4 .* RF ./ (2 * BGrads .* Constants.zeemansplit),'--','Color', [0.4 0.4 0.4]);
xlabel('Quadrupole gradient (G/cm)', 'Interpreter', 'Latex');
ylabel('Atom position ($\mu$m)', 'Interpreter', 'Latex');
title('Position below quadrupole');
legend([actualH theoryH], 'actual', 'resonant ellipsoid', 'Location', 'SouthEast');

%%

%%
% plot(BTOPs,tms); hold on;
% dz = ((resonantEllipsoidWidth(1.4, BGrad, zsf).^2 - (BTOPs/(BGrad*1e-4)).^2).^0.5)/2;
% plot(BTOPs, -dz, ':'); hold off
