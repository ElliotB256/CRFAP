function [ F ] = sortByCurvature( x, F )
%SORTBYCURVATURE Sorts the given energy levels by their curvature.
% Sorts energy levels by taking the second derivative. 'x' is used as a
% generalised x coordinate, which can be spatial, zeeman splitting etc.
%
% SYNTAX: sortByCurvature( x, F )

dx = x(2:end) - x(1:end-1);
dydx  = (F(:,2:end) - F(:, 1:end-1))./repmat(dx, size(F,1), 1);
dx = (dx(1:end-1) + dx(2:end))/2;
d2ydx2 = (dydx(:,2:end) - dydx(:,1:end-1))./repmat(dx, size(dydx,1), 1);

% Sum curvature
s = sum(d2ydx2, 2);

% order energy levels
[~,j] = sort(s);

F = F(j,:);


end

