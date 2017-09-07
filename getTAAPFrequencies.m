function [ result ] = getTAAPFrequencies( RF, BTOP, RFAmp, QuadGrad, varargin )
%GETTAAPFREQUENCIES Calculates trap frequencies for a TAAP trap.
%   Assumes the RF is circularly polarised.
%
%   RF: Dressing RF frequency (MHz)
%
%   BTOP: TOP field amplitude, Gauss
%
%   RFAmp: RF Amplitude, B, in Gauss.
%
%   QuadGrad: Quadrupole gradient in Gauss/cm. Defined such that B_z = 2 *
%   B' z.
%
%   Optional input arguments:
%    'Species': '85' or '87', lower hyperfine states.
% 
%   Syntax: getTAAPFrequencies( RF, BTOP, RFAmp, QuadGrad, ... )

import Constants.*

ip = inputParser();
ip.addParameter('Species', 87);
ip.parse(varargin{:});

% Species: 85 or 87
switch ip.Results.Species
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

% Get displacement of quadrupole centre in microns
xTOP = BTOP/QuadGrad * 1e4;

% Define potential
trap = @(a,b,c,t) SRF.ShellTrap(...
    a - xTOP .* cos(2*pi*t), b - xTOP .* sin(2*pi*t),c,...
    zsf, QuadGrad, RF, RFAmp, mFtilde) + ...
    gpe(c, mass);

% Perform time averaging crudely along line x,y=0 to find trap
% minimum.
z = 0:-.1:-2000;
x = zeros(size(z));
y = x;

pot = timeAverage(@(t) trap(x,y,z,t), 10);

[~,ip] = min(pot);
trapMinZ = z(ip);

clear ip pot;

% Time average along x,y,z axis at the minimum spot.
% Take a small size around the trap center.
probeLength = 2; %um
probeResolution = 100;
ps = -probeLength:probeLength*2/probeResolution:probeLength;
zs = zeros(size(ps));
Nt = 10;

trapFx = Util.getTrapFreq(ps, timeAverage(@(t) trap(ps,zs,zs+trapMinZ,t),Nt), mass);
trapFy = Util.getTrapFreq(ps, timeAverage(@(t) trap(zs,ps,zs+trapMinZ,t),Nt), mass);
trapFz = Util.getTrapFreq(ps+trapMinZ, timeAverage(@(t) trap(zs,zs,ps+trapMinZ,t),Nt), mass);


result = struct('fx', trapFx(2), 'fy', trapFy(2), 'fz', trapFz(2));

end

