function [ minima ] = getMinima( x, y )
%GETMINIMA Surveys a given level y and extracts local minima. Returns
%minima x.
% SYNTAX: getMinima(x, y)

dydx = derivative(y, x);
d2ydx = derivative(dydx, x);

% Locate points of interest where the derivative intersects zero. Identify
% as minima those where the second derivative is positive.

i = find(abs(diff(sign(dydx))) > 0);
% plot(x,y); hold on; plot(x(i), y(i),'r.'); hold off;

% intersect at zero occurs between ith and i+1th element
isct = [i' i'+1];

avgQty = @(b,a) (b(a(:,1)) + b(a(:,2)))/2;
minj = i(avgQty(d2ydx, isct) > 0);

minima = avgQty(x, [minj' minj'+1]);
%plot(x,y); hold on; plot(minima*[1 1], ylim,'r-'); hold off;

%TODO: iterative search to find minimum of x using interpolated function or
%fitted parabola. Probably not neccessary right now - can always turn up
%the number of interpolations on the MRF generation to get the same effect.

end

