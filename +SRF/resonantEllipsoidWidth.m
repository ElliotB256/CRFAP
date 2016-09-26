function [ rew ] = resonantEllipsoidWidth( RF, BGrad, zsf )
%RESONANTELLIPSOIDWIDTH Calculates the width of the resonant ellipsoid
% RF: Radiofrequency in MHz
% BGrad: Quadrupole field gradient in Gauss/cm
% zsf: Zeeman splitting factor in MHz/Gauss

rew = RF / zsf / BGrad * 1e4;


end

