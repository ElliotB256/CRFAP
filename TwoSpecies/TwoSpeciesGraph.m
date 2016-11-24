%% Two Species
% Create a graph showing the potentials we could make for Rb 85 and 87.
% This will trap the 87 in the potential due to 3.5 MHz, and 85 in the
% potential due to 4.0/.1/.2.

n = [ 35 36 37 38 39 40 41 42 43 44 45 ]'; % prev was just to 42
df = 0.1;
RFs = n*df;
Rabi = df * [ 0.7 0 0 0 0 0 0 0 0.7 0.3 0.7 ]';
qdrpGrad = 200;
Bs = 3.4:0.05:4.7;

df = 0.1;
% 85 single well, 87 double well
% n =         [ 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 ]';
% Rabi = df * [ .7 00 00 00 00 00 00 00 00 00 00 00 00 00 .7 .4 .7 ]';
n =         [ 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 ]';
Rabi = df * [ .7 .4 .7 00 00 00 00 00 00 00 00 00 00 00 00 00 .7 ]';
RFs = n*df;
Bs = 2.7:0.05:4.7;

warning('This may take a while...');
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 6, 'qdrpGrad', qdrpGrad);
%F = MRF.sortEnergies(B, MRF.ladder(RFs, 10, F));
%lad = MRF.ladder(RFs, 3, F);
%glad = lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1);

plot(B, F);


%%
% A very large number of avoided crossings result. Does that look to cause
% a problem for us? Try performing the same calculation using just two
% dressing rfs. For comparison purposes, I've decided to keep the df the
% same - this is in some way a test of how accurate the numerics are.

n = [ 35 36 37 38 39 40 41 42 ]';
df = 0.1;
RFs = n*df;
Rabi = [ 0.2 0 0 0 0 0.0 0.0 0.2 ]';
qdrpGrad = 200;

Bs = 3.4:0.05:4.4;
warning('This may take a while...');
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 6, 'qdrpGrad', qdrpGrad);
F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 10, F));
lad = MRF.ladder(RFs, 3, F2);
glad = lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1);

plot(B, glad);

%%
% Plot section 1 potential in a clearer way.

figure(1);
subplot(1, 4, 1:3);
lad = MRF.ladder(RFs, 10, F);
plot(B, (lad-0.5)*1e3, 'k');
ylim([-0.01 .11 ] * 1e3);
xlim([2.8 4.7]);
xlabel('Zeeman splitting (MHz)', 'Interpreter', 'Latex');
ylabel('Energy (kHz)', 'Interpreter', 'Latex');

% Add grey circles to mark finite depths of wells
hold on;
plot(3.163, 89.54, 'o', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 1.5, 'MarkerSize', 11);
plot(2.838, 89.57, 'o', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 1.5, 'MarkerSize', 11);
plot(4.424, 98.66, 'o', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 1.5, 'MarkerSize', 11);
plot(4.568, 98.66, 'o', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 1.5, 'MarkerSize', 11);

% Add dashed lines to indicate the dressing rf components
plot([ 1 1 ] * 2.9, ylim, ':', 'Color', [ 1 1 1 ] * 0.7, 'LineWidth', 1.5);
plot([ 1 1 ] * 3.0, ylim, ':', 'Color', [ 1 1 1 ] * 0.7, 'LineWidth', 1.5);
plot([ 1 1 ] * 3.1, ylim, ':', 'Color', [ 1 1 1 ] * 0.7, 'LineWidth', 1.5)
plot([ 1 1 ] * 4.5, ylim, ':', 'Color', [ 1 1 1 ] * 0.7, 'LineWidth', 1.5);

hold off;
ax1 = gca;
set(gca,'ActivePositionProperty','outerposition');

subplot(1, 4, 4);
% plot trapped manifolds for 85 and 87
trap = F(1,:)+0.1;

gF87 = 0.7;
gF85 = 2/3 * gF87;

plot(B / gF87 * (1e4/ qdrpGrad) - 322 + 0.404, trap * 1e3, 'Color', [ 0.8 0.4 0.3], 'LineWidth', 1.5); hold on;
plot(B / gF85 * (1e4/ qdrpGrad) - 322 + 0.404, trap * 1e3, 'Color', [ 0.3 0.3 0.8], 'LineWidth', 1.5); hold off;
xlim( [ -12 12 ]);
%legend({'$^{87}$Rb', '$^{85}$Rb'}, 'Interpreter', 'Latex');
set(gca, 'YTick', 50:10:100);

% Link axes together.
% set(gca, 'Units', 'Pixels');
% set(ax1, 'Units', 'Pixels');

%pos = get(gca, 'Position');
ylim([ 50 100 ]);

xlabel('z ($\mu$m)', 'Interpreter', 'Latex');
ylabel('Energy (kHz)', 'Interpreter', 'Latex');
ax2 = gca;
set(gca,'ActivePositionProperty','outerposition');

% Final figure adjustments
set(gcf, 'Color', [ 1 1 1 ]);
set(gcf, 'Units', 'pixels', 'Position', [ 200 200 800 300 ]);
set(ax1, 'Units', 'pixels');
set(ax2, 'Units', 'pixels');
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos1 = [ pos1(1)-50 pos1(2) pos1(3)+50 pos1(4) ];
pos2 = [ pos2(1)+50 pos2(2) pos2(3) pos2(4) ];
set(ax1, 'Position', pos1);
set(ax2, 'Position', pos2);
set(ax1, 'Units', 'Normalized');
set(ax2, 'Units', 'Normalized');

set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');

w = pos(3);
h = pos(4);
p = 0.01;
set(gcf,...
  'PaperUnits','centimeters',...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);
