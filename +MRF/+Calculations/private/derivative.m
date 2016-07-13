function [ dydx ] = derivative( y, x )
%DERIVATIVE Differentiates y with respect to x
% Returns a vector of length == length(y).

dy = diff(y);
dy = ([dy dy(end)] + [dy(1) dy])/2;

dx = diff(x);
dx = ([dx dx(end)] + [dx(1) dx])/2;

dydx = dy ./ dx;

end

