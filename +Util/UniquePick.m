function [ b ] = UniquePick( a, b, thresh )
%UNIQUEDBL Returns all elements from b that are not in a and only each
%element once.

if nargin < 3
thresh = min(diff(sort(unique(a))))/10;
end

mask = ones(size(b));
for i=1:length(b)
   d = abs(repmat(b(i), size(a, 1), size(a, 2)) - a);
   if any(d < thresh)
       mask(i) = 0;
   end
end
b(~mask) = [];

% Perform unique sweep over resulting b
b = uniquetol(b, thresh)';

end