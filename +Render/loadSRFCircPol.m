function [ theta, x, z, v, sim ] = loadSRFCircPol( filename )
%LOADSRFCIRCPOL Loads a file generated by SRFCircPol
% Returns 4 vectors corresponding to value of theta, x, z, potential for
% each point in space that data exists for. Also returns a struct, sim,
% describing the parameters of the trap.

warning('deprecated: use loadCircPol');
[ theta, x, z, v, sim ] = Render.loadCircPol( filename );

end

