function [ gp ] = gpe( ZeemanSplit, qdrpGrad )
%GPE Calculates gravitational potential energy in MHz. Assumes the field
%points are on the vertical axis. The bare zeeman splitting (ZeemanSplit)
%is mapped to spatial coordinates using the given quadrupole gradient.
% Syntax: gpe( ZeemanSplit, qdrpGrad )
%  ZeemanSplit: energy splitting of the bare zeeman states in the
%               quadrupole
%  qdrpGrad: quadrupole gradient in G/cm

BGauss = ZeemanSplit / Constants.gF;
microns = -BGauss / qdrpGrad * 1e4;
gp = gpe(microns, 87);

end

