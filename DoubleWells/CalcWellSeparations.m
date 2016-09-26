% Calculate well separations over a range of different frequency
% separations and quadrupole gradients.

Bs = RFs(1):(RFs(3)-RFs(2))/3:1.1*(RFs(3)-RFs(1))+RFs(1);

[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 4, 'qdrpGrad', qdrpGrad);
F = MRF.sortEnergies(B, MRF.ladder(RFs, 3, F));
lad = MRF.ladder(RFs, 3, F);
glad = lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1);

% Find the well separation for the lowest energy level.
Fi = @(x) interp1(B, glad(1, :), x);
x1 = fminbnd(Fi, RFs(1), RFs(2));
x2 = fminbnd(Fi, RFs(2), Bs(end));

plot(B, glad(1,:), '.-', 'Color', [0.3 0.3 0.3]);
hold on; plot(B, Fi(B), ':', 'Color', [0.0 0.0 0.0]); hold off
hold on; plot(x1, Fi(x1), '.', 'MarkerSize', 14, 'Color', [0.8 0.3 0.3]); hold off;
hold on; plot(x2, Fi(x2), '.', 'MarkerSize', 14, 'Color', [0.8 0.3 0.3]); hold off;

% Determine the distance between the wells.
wellSep = (x2 - x1)/qdrpGrad ...% in cm
            * 1e4; % in um
fprintf('Wells with qdrpGrad=%.0f G/cm and Delta_f=%0.2f have spatial separation of %.1f um.\n', qdrpGrad, df, wellSep)

% Illustrate separation on the digram
hold on;
my = mean([Fi(x1), Fi(x2)]);
yl = ylim;
plot([x1, x2], [my my], '-', 'Color', [0.5 0.5 0.5]);
plot([x1 x1], yl, ':', 'Color', [0.5 0.5 0.5]);
plot([x2 x2], yl, ':', 'Color', [0.5 0.5 0.5]);
t = text(mean([x1 x2]), my, [num2str(wellSep, '%.1f') '$\mu m$'], 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
hold off;
ylim(yl);

title(['$B''= ' num2str(qdrpGrad, '%.1f') '$ G/cm, $\Delta_f=' num2str(df, '%.2f') '$ MHz'], 'Interpreter', 'Latex');

% Calculate the trap frequencies for both wells