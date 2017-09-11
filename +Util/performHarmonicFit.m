function [ k2 ] = performHarmonicFit( x, y )
%PERFORMHARMONICFIT Performs a fit and returns quadratic parameters within 95%
%confidence interval


%Rescale and center the data points
xN = (x - mean(x))/std(x);
yN = (y - mean(y))/std(y);

f = fit(xN', yN', 'poly2');
ci = confint(f);

k2 = std(y)/std(x).^2 .* [ci(1,1) f.p1 ci(2,1)];

%Quick test, note: Won't remove gradient (~x terms)
% figure(2)
%    [~,i] = min(y);
%    xt = x - x(i);
%    yt = k2(2) .* xt.^2;
%    plot(x-x(i),y-y(i),'x'); hold on; plot(xt,yt); hold off
% figure(1)

end