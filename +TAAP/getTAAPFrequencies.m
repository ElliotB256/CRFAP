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

TC = TAAP.Calculator(RF, RFAmp, BTOP, QuadGrad);
TC.WithSpecies(ip.Results.Species, 'mFTilde', ip.Results.mFTilde);
TC.WithAnisotropy(ip.Results.Anisotropy);

% Perform time averaging crudely along line x,y=0 to find trap
% minimum.
z = 0:-.1:-2000;
x = zeros(size(z));
y = x;

% pot = Util.timeAverage(@(t) trap(x,y,z,t), ip.Results.TimeAverageSteps);
pot = TC.Calculate(x,y,z);
[~,ip] = min(pot);

if (ip == 1 || ip == length(z))
    plot(z, pot);
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

trapFx = Util.getTrapFreq(ps, TC.Calculate(ps,zs,zs+trapMinZ), TC.Atom.Mass);
trapFy = Util.getTrapFreq(ps, TC.Calculate(zs,ps,zs+trapMinZ), TC.Atom.Mass);
trapFz = Util.getTrapFreq(ps+trapMinZ, TC.Calculate(zs,zs,ps+trapMinZ), TC.Atom.Mass);

result = struct('fx', trapFx(2), 'fy', trapFy(2), 'fz', trapFz(2), 'sag', trapMinZ);

end