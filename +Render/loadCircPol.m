function [ theta, x, z, v, sim ] = loadCircPol( filename )
%LOADSRFCIRCPOL Loads a file generated by CircPol
% Returns 4 vectors corresponding to value of theta, x, z, potential for
% each point in space that data exists for. Also returns a struct, sim,
% describing the parameters of the trap.

fH = fopen(filename, 'r');
oc = onCleanup(@() fclose(fH));

n = fread(fH, 1, 'uint32');
rf = fread(fH, n, 'double');
n = fread(fH, 1, 'uint32');
rabi = fread(fH, n, 'double');
qdrpGrad = fread(fH, 1, 'double');

sim.rf = rf; sim.rabi = rabi; sim.qdrpGrad = qdrpGrad;

    function [th, x,z,v] = readStrip(fH)
        th = fread(fH, 1, 'double');
        if ~isempty(th);
            N = fread(fH, 1, 'double');
            x = fread(fH, N, 'double');
            z = fread(fH, N, 'double');
            nV = fread(fH, 1, 'uint32');
            for i=1:nV
                v(:,i) = fread(fH, N, 'double');
            end
        else
            N = []; x = []; z = []; v = [];
        end
        
        th = repmat(th, size(x));
    end

theta = []; x = []; z = []; v = [];
while ~feof(fH)
    [tth, tx, tz, tv] = readStrip(fH);
    theta = [ theta tth ];
    x = [ x tx ];
    z = [ z tz ];
    v = cat(3, v, tv);
end

% remove non singleton dimensions for ease of use when only one level is required.
v = squeeze(v);

end

