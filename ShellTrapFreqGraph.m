%% Produces a nice graph of shell trap frequencies
% Start by precalculating required information...

% Calculate shell trap frequencies for a range of quadrupole gradients and
% dressing RF amplitudes.

RF = 4.2; %MHz
BGrads = 10:5:600;
RFAmps = (0.1:0.05:0.7)./Constants.zeemansplit;

freqs = zeros(3,length(RFAmps),length(BGrads));

%for loop yuck
for i=1:length(BGrads)
    for j=1:length(RFAmps)
        
        f = shellTrapFrequencies(RF, RFAmps(j), BGrads(i));
        freqs(:,j,i) = f;
        
    end
end

%% Create the graph itself

figure; ax1 = axes(); hold on

% plot axial trap frequencies in shades of blue. Each series corresponds to
% different BRF.
for i=1:length(RFAmps)
    lA = RFAmps(i)/max(RFAmps(:));
    palette = [0.2 0.2 .6; .6 .6 1];
    color = palette(1,:) .* lA + (1-lA) .* palette(2,:);
    ah = plot(BGrads, squeeze(freqs(3, i, :)), 'Color', color, 'LineWidth', 1+lA*0.5); 
end

set(ax1, 'YAxisLocation', 'left')

% create second axis object
 ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none'); hold on
set(ax2, 'XLim', xlim(ax1));

% plot radial trap frequencies in shades of red.
for i=1:length(RFAmps)
    lA = RFAmps(i)/max(RFAmps(:));
    % need to mask it so we don't plot antitrapping regions
    mask = imag(squeeze(freqs(1,i,:))) < eps;
    ind = find(~mask, 1, 'first');
    if ~isempty(ind)
        mask(ind) = 1;
        freqs(1,i,ind) = 0;
    end
    clear ind;
    
    radialF = squeeze(freqs(1, i, :));
    palette = [.8 0.2 .2; 1 0.8 0.8];
    color = palette(1,:) .* lA + (1-lA) .* palette(2,:);
    rh = plot(BGrads(mask), radialF(mask), 'Color', color, 'LineWidth', 1+lA*0.5); 
end

% Label axes
xlabel(ax1, 'Quadrupole Gradient (G/cm)', 'Interpreter', 'Latex');
ylabel(ax1, 'Axial (Hz)', 'Interpreter', 'Latex');
ylabel(ax2, 'Radial (Hz)', 'Interpreter', 'Latex');
legend([ah rh], 'Axial', 'Radial', 'Location', 'NorthWest');

title(['Single-RF Shell Trap (' num2str(RF, '%.1f') ' MHz)']);
tidyFigure(0.1);
saveas(gcf, 'ShellTrapFreqs.pdf');


%% Another useful figure: Ratio of freqs
% This figure shows the ratio of trap freqs, aka anisotropy of trap.
figure; hold on
for i=1:length(RFAmps)
    lA = RFAmps(i)/max(RFAmps(:));
    %palette = [.8 0.2 .4; 0.2 0.2 0.8];
    palette = [.8 0.2 .2; 1 0.8 0.8];
    color = palette(1,:) .* lA + (1-lA) .* palette(2,:);
    mask = abs(imag(squeeze(freqs(1,i,:)))) < eps;
    ratio = squeeze(freqs(3,i,:))./squeeze(freqs(1,i,:));
    plot(BGrads(mask), ratio(mask), 'Color', color, 'LineWidth', 1+lA*0.5);
end

xlabel('Quadrupole Gradient (G/cm)', 'Interpreter', 'Latex');
ylabel('$\omega_z / \omega_x$', 'Interpreter', 'Latex', 'FontSize', 16);
title(['Shell Trap Anisotropy (' num2str(RF, '%.1f') ' MHz)']);
tidyFigure(0.1);
saveas(gcf, 'ShellTrapFreqRatio.pdf');
