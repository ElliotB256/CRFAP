function [ f ] = shellTrapFrequencies( RF, BRF, BGrad )
% SHELLTRAPFREQUENCIES Calculates trap frequencies for a single RF dressed
% AP.
%  BRF: Amplitude of dressing RF in MHz
%  RF: RF frequency in MHz
%  BGrad: Qdrp gradient, Gauss/cm
import Constants.*

zsf = Constants.zeemansplit; %Mhz/Gauss
rb = rubidium87();
mass = rb.mass; %amu

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

clear ip pot;

% Time average along x,y,z axis at the minimum spot.
% Take a small size around the trap center.
probeLength = 5; %um
probeResolution = 100;
ps = -probeLength:probeLength*2/probeResolution:probeLength;
zs = zeros(size(ps));

f = zeros(1,3);

fx = getTrapFreq(ps,trap(ps,zs,zs+trapMinZ), mass);
fy = getTrapFreq(ps,trap(ps,zs,zs+trapMinZ), mass);
fz = getTrapFreq(ps+trapMinZ,trap(zs,zs,ps+trapMinZ), mass);

f = [fx(2) fy(2) fz(2)];

end
