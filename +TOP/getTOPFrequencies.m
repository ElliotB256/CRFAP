function [ f ] = getTOPFrequencies( B, Bgrad, mu, m )
%GETTOPFREQUENCIES Calculate TOP trap frequencies
% 
%  B is a vector defining the magnitude of the TOP field in Gauss. If
%  length(B) == 1, it is taken that B_x = B_y, otherwise we take B_x = B(1)
%  and B_y = B(2).
%  
%  Bgrad is the quadrupole gradient in the radial plane, Gauss/cm. B =
%  Bgrad(x,y,-2z).
%  
%  mu is the magnetic dipole of the atom. It is equal to mFgFmuB. The
%  number here should be specified in MHz/G.
%  
%  m is the mass in atomic mass units
%  
%  Syntax: f = getTOPFrequencies(B, Bgrad, mu, m)

if length(B) == 1
   B = [ B B ];
end

% Create a function to calculate the instantaneous potential energy.
% Answer is returned in MHz for mu in MHz/G
U = @(x,y,z,nt) mu * Bgrad * ( ...
    (1e-4.*x + B(1)/Bgrad .* cos(2*pi*nt)).^2 ...
    + (1e-4.*y + B(2)/Bgrad .* sin(2*pi*nt)).^2 ...
    + (-2 .* 1e-4.*z).^2).^0.5;

% Determine lengths for probing the potential over. Lengths are in um.
p = linspace(-10, 10, 30);
a = zeros(size(p));

% Time average the potential over one period
calcU = @(x,y,z) Util.timeAverage(@(nt) U(x,y,z,nt), 100);
potx = calcU(p,a,a);
poty = calcU(a,p,a);
potz = calcU(a,a,p);

fx = Util.getTrapFreq(p, potx, m);
fy = Util.getTrapFreq(p, poty, m);
fz = Util.getTrapFreq(p, potz, m);
f = struct('fx', fx(2), 'fy', fy(2), 'fz', fz(2));

end

