function [ newX ] = Refine( x, y )
%REFINE Calculates refinement points for the given mesh

if (size(y,1)) ~= 1
    y = y';
end

if (size(y,1)) ~= 1
    error('only works for vectors');
end

if (size(x,1)) ~= 1
    x = x';
end

if (size(x,1)) ~= 1
    error('only works for vectors');
end

% identify any points where the derivative changes sign
d = y - [y(1) y(1:end-1)];
s = sign(d);

% stationary point somewhere between j1 and j2
j1 = find(abs(diff(s)) > 0.1);
if length(j1) > 1
    j1 = j1(2:end);
    j2 = min(j1+1,length(d));
    newX = mean([x(j1); x(j2)], 1);
    newX = [newX mean([x(j1-1); x(j2-1)], 1)];
    newX = newX(:);
else
    newX=[];
end

newX = unique(newX);

% no elements of newX should exist in current set
% thresh = min(diff(sort(unique(x))))/10;
% mask = ones(size(newX));
% for i=1:length(newX)
%    d = abs(repmat(newX(i), size(x, 1), size(x, 2)) - x);
%    if any(d < thresh)
%        mask(i) = 0;
%    end
% end
% newX(~mask) = [];

newX = Util.UniquePick(x, newX);

end

