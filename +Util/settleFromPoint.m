function [ minloc ] = settleFromPoint( x, y, x0 )
%SETTLEFROMPOINT Finds a minimum by settling from a given start point into
%a local minimum.
% SYNTAX: settleFromPoint(x, y, x0)

% Find position of x0 along x
minloc = x0;

i = find(x0 < x, 1, 'first');

% x0 is part way between i-1 and i
if isempty(i)
    error('x0 not within range x');
end

settled = 0;
previ = NaN;

while ~settled
    
    if i == 1 || i == length(y)
        if i==1 && y(i) > y(i+1)
            i = i + 1;
        elseif i==length(y) && y(i) > y(i-1)
            i = i - 1;
        else
            error('Reached edge of potential during minimisation.');
        end
    end
    
    if y(i-1) < y(i)
        i = i - 1;
    elseif y(i+1) < y(i)
        i = i + 1;
    end
    
    if previ == i
        settled = 1;
    end
    
    previ = i;
end

minloc = i;

end

