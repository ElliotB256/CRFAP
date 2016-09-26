function [ gp ] = gpe( ZeemanSplit, Bprime, varargin )
%GPE Calculates gravitational potential energy in MHz. Assumes the field
%points are on the vertical axis. The bare zeeman splitting (ZeemanSplit)
%is mapped to spatial coordinates using the given quadrupole gradient. Note
%that we use the convention that B_z = 2 * B' * z!
% Syntax: gpe( ZeemanSplit, qdrpGrad )
%  ZeemanSplit: energy splitting of the bare zeeman states in the
%               quadrupole
%  qdrpGrad: quadrupole gradient in G/cm

p = inputParser;
   addRequired(p,'ZeemanSplit',@isnumeric);
   addRequired(p,'Bprime',@isnumeric);
   addParameter(p,'gF', Constants.gF, @isnumeric);
   addParameter(p,'mass', 87, @isnumeric);
   
parse(p, ZeemanSplit, Bprime, varargin{:});
   
BGauss = ZeemanSplit / p.Results.gF;
microns = -BGauss / Bprime * 1e4 / 2;

% above calculation verified as correct height.
gp = gpe(microns, p.Results.mass);

end