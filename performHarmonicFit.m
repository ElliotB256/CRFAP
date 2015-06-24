function [ k2 ] = performHarmonicFit( x, y )
%PERFORMHARMONICFIT Performs a fit and returns quadratic parameters within 95%
%confidence interval

f = fit(x', y', 'poly2');
ci = confint(f);

k2 = [ci(1,1) f.p1 ci(2,1)];
end

