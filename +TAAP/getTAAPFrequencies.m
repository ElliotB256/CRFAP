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
ip.addParameter('mFTilde', []);
ip.addParameter('TimeAverageSteps', 20);
ip.addParameter('Anisotropy', 1); % B_TOP y = anisotropy * B_TOP x
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

if ~isempty(ip.Results.mFTilde)
   mFtilde = ip.Results.mFTilde; 
end

% Get displacement of quadrupole centre in microns
xTOP = BTOP/QuadGrad * 1e4;

% Define potential
trap = @(a,b,c,nt) SRF.ShellTrap(...
    a - xTOP .* cos(2*pi*nt), ...
    b - ip.Results.Anisotropy * xTOP .* sin(2*pi*nt), ...
    c,...
    zsf, QuadGrad, RF, RFAmp, mFtilde) + ...
    Util.gpe(c, mass);

% Perform time averaging crudely along line x,y=0 to find trap
% minimum.
z = 0:-.1:-2000;
x = zeros(size(z));
y = x;

pot = Util.timeAverage(@(t) trap(x,y,z,t), ip.Results.TimeAverageSteps);

[~,ip] = min(pot);

if (ip == 1 || ip == length(z))
    error('Cannot find trap minimum in vertical direction; reached end of range');
end
trapMinZ = z(ip);

clear ip pot;

% Time average along x,y,z axis at the minimum spot.
% Take a small size around the trap center.
probeLength = 2; %um
probeResolution = 100;
ps = linspace(-probeLength, probeLength, probeResolution);
zs = zeros(size(ps));
Nt = 10;

trapFx = Util.getTrapFreq(ps, Util.timeAverage(@(t) trap(ps,zs,zs+trapMinZ,t),Nt), mass);
trapFy = Util.getTrapFreq(ps, Util.timeAverage(@(t) trap(zs,ps,zs+trapMinZ,t),Nt), mass);
trapFz = Util.getTrapFreq(ps+trapMinZ, Util.timeAverage(@(t) trap(zs,zs,ps+trapMinZ,t),Nt), mass);


result = struct('fx', trapFx(2), 'fy', trapFy(2), 'fz', trapFz(2), 'sag', trapMinZ);

end