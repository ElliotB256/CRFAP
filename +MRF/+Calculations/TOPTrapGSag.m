%% TOP Trap Gravitational Sag
% Calculating the gravitational sag in the TOP trap used for measurements
% of the quadrupole position.

qdrpI = 200; %200 A Qdrp at end of standard TOP bec sequence.
qdrp  = 62.5 / 20 * qdrpI; % Gauss/cm

TOPControlV = 5; %V on experimental control VI
%BTopI  = 7.79 * TOPControlV; % Take TOP calibration from 2016_05_11

% temporary calibration of TOP: 1V on vi = 6 MHz splitting for imaging
TOP = TOPControlV * 6 / 0.7; % Gauss

% qdrp potential on z axis. Due to cylindrical symmetry of TOP the z-value
% won't change.
z = -200:1:200; % microns;
x = zeros(size(z)); y = x;

field = qdrp * [-x/2; -y/2; z];