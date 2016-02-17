function tidyFigure(padding)
%TIDYFIGURE Tidies a figure to save without whitespace

if nargin < 1
    padding = 0.6;
end

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');
li = get(gca,'LooseInset');
ti = li*padding + ti.*(1-padding);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

end