% Simulate single frequency shell
clear;
import Constants.*

zsf = Constants.zeemansplit; %Mhz/Gauss
BGrad = 84; %Gauss/cm
rb = rubidium87();
mass = rb.mass; %amu
BRF = 0.5;
RF = 4;

tms = [];

BTOP = 0;
   
BGrads = 30:5:300;

for BGrad=BGrads
    
    % Define potential
    trap = @(a,b,c) ShellTrap(...
        a, b, c,...
        zsf, BGrad, RF, BRF) + ...
        gpe(c, mass);
    
    % first iter: crudely find trap minimum.
    z = 0:-1:-1000;
    x = zeros(size(z)); y = x;
    pot = trap(x,y,z);
    
    [~,ip] = min(pot);
    trapMinZ = z(ip);

    % second iter: refine trap minimum
    z = (-3:0.001:3) + trapMinZ;
    x = zeros(size(z)); y = x;
    pot = trap(x,y,z);
    [~,ip] = min(pot);
    trapMinZ = z(ip);
    tms(end+1) = trapMinZ;
    
    clear ip pot;
end

%%

figure; ax1 = axes();
actualH = plot(BGrads, tms, 'LineWidth', 1, 'Color', [0.4 0.4 1]); hold on
theoryH = plot(BGrads, -1e4 .* RF ./ (2 * BGrads .* Constants.zeemansplit),'--','LineWidth', 1,'Color', [0.4 0.4 0.4]);
xlabel('Quadrupole gradient (G/cm)', 'Interpreter', 'Latex');
ylabel('Atom position ($\mu$m)', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', [0.4 0.4 1]);
title({'Position below quadrupole', ['RF = ' num2str(RF, '%.1f') 'MHz, Amplitude = ' num2str(BRF, '%.1f') 'MHz']});

% create second axis object
 ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none'); hold on
set(ax2, 'XLim', xlim(ax1));
sag = plot(BGrads, tms + 1e4 .* RF ./ (2 * BGrads .* Constants.zeemansplit), '--', 'LineWidth', 1.5, 'Color', [1 0.4 0.4]);
legend([actualH theoryH sag], 'Actual Position', 'Resonant Ellipsoid Position', 'Sag','Location', 'SouthEast');
ylabel('Sag ($\mu$m)', 'Interpreter', 'Latex', 'Color', [.9 0.4 0.4], 'FontSize', 14);