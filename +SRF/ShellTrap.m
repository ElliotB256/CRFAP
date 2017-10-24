function [ pot ] = ShellTrap( x,y,z, zsf, BGrad, RF, RFAmp, mFtilde, BRF )
%SHELLTRAP 

% x,y,z : coordinates of field points, in microns
%  
% zsf   : Zeeman split factor, g_f \mu_B / hbar, in units of MHz/Gauss.
% 
% BGrad : Quadrupole gradient, Gauss/cm. Defined such that B_z = 2 B' z
% 
% BRF   : a function handle that returns the amplitude of the perpendicular
%         component of  applied RF. Defaults to circ polarisation along z,
%         with node at the top.
% 
% RFAmp : Amplitude of dressing RF in Gauss
% 
% Returns: potential energy at field points in MHz

if nargin < 8
    mFtilde = 1;
end

%Default to circular polarisation
if nargin < 9
   BRF = @(a,b,c) ( 1 - 4 .* c ./ (a.^2 + b.^2 + 4 .* c.^2).^0.5 + ...
       4 .* c.^2 ./ (a.^2 + b.^2 + 4 .* c .^2)).^0.5;
end

%Calculate magnitude of B field in Gauss.
Bmod = BGrad.*(x.^2 + y.^2 + 4.*z.^2).^0.5 .* 1e-4;

pot = mFtilde * ((zsf .* Bmod - RF).^2 + (zsf/2 .* RFAmp .* BRF(x,y,z)).^2).^0.5;

end

