function [ labels ] = labelStates( N )
%LABELSTATES Labels states in the bare basis
% N(0) should be number of bare atomic levels.
% N(1+) should be max number of photons in mode.

labels = int32(zeros(prod(N(:)), length(N)));

for i=1:prod(N(:))
    [temp{1:length(N)}] = ind2sub(N, i);
    labels(i,:) = int32([temp{:}]);
end

end

