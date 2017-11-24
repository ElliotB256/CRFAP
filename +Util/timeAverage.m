function [ tap ] = timeAverage( potential, N )
%TIMEAVERAGE Time averages the given potential with the specified number of
%steps.
% potential: the potential to average N+1 times. Should be a function handle
% which accepts one argument, normalised time in the range 0-1.

if N < 1
    error('N must be a positive integer');
end

nt = 0:1/N:1;
% remove final point - otherwise we sample that point twice!
nt(end) = [];

for t=nt
    ip = potential(t);
    if ~exist('tap', 'var')
        tap = zeros(size(ip));
    end
    tap = tap + ip;
end

tap = tap/(N+1);

end