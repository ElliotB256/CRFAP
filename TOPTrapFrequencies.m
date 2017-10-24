function [ s ] = TOPTrapFrequencies( Btop, BGrad, mF, uBgF, m )
%TOPTRAPFREQUENCIES Gets frequencies of the TOP trap
% Btop: TOP Field, Gauss
% BGrad: Quadrupole gradient, G/cm, such that B_z = 2 BGrad z
% m: mass, amu
% mF: magnetically trapped state
% uBgF: Lande g-factor times Bohr magneton, MHz/Gauss

%Pethick and Smith, 4.1.2 p63
fx = (Constants.hbar * 2 * pi * uBgF * 1e6 * abs(mF) * (BGrad).^2 / (2 * m * Constants.amu * Btop) * 1e4).^0.5 / 2 / pi;
fz = sqrt(8) * fx;
s = struct('fx', fx, 'fy', fx, 'fz', fz, 'wx', 2*pi*fx, 'wy', 2*pi*fx, 'wz', 2*pi*fz);

end

