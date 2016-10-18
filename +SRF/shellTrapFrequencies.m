function [ f ] = shellTrapFrequencies( RF, BRF, BGrad )
% SHELLTRAPFREQUENCIES Calculates trap frequencies for a single RF dressed
% AP.
%  BRF: Amplitude of dressing RF in Gauss
%  RF: RF frequency in MHz
%  BGrad: Qdrp gradient, Gauss/cm
% 
% Note: for very shallow BGrad the bottom of the shell will extend beyond
% the meshing region and we will have erroneous results!
import Constants.*
import SRF.*

zsf = Constants.zeemansplit; %Mhz/Gauss
rb = rubidium87();
mass = rb.mass; %amu

% Define potential
trap = @(a,b,c) ShellTrap(...
    a, b, c,...
    zsf, BGrad, RF, BRF) + ...
    gpe(c, mass);

% first iter: crudely find trap minimum.
z = 0:-10:-10000;
x = zeros(size(z)); y = x;
pot = trap(x,y,z);

[~,ip] = min(pot);
trapMinZ = z(ip);

% second iter: refine trap minimum
z = (-30:0.001:20) + trapMinZ;
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

fx = Util.getTrapFreq(ps,trap(ps,zs,zs+trapMinZ), mass);
fy = Util.getTrapFreq(ps,trap(ps,zs,zs+trapMinZ), mass);
fz = Util.getTrapFreq(ps+trapMinZ,trap(zs,zs,ps+trapMinZ), mass);

f = [fx(2) fy(2) fz(2)];

end

