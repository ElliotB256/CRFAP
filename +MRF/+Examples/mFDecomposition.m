%% mF Decomposition
% Calculate the mF composition of eigenvectors of the propagator T, the
% application of which evolves the wavefunction for one period of time. As
% this operator T commutes with the Hamiltonian, these eigenvectors are
% also the eigenvectors of the dressed Hamiltonian.
RF = 3;
Rabi = 0.3;
qdrpGrad = 100;

Bs=linspace(2.7, 4.5, 40);
[ eigF, eigV ] = MRF.GetQuasiEnergies(Bs, RF, Rabi, 'parallel', 0);

% Having calculated the mF amplitudes, lets convert these into a color map.
% We map the abs() of each mF state amplitude into the colors r,g,b. The
% result is plotted using a patch object so that we can vary color along
% the line. Each state is plotted in turn.
for i=1:3
    x = [ Bs NaN ]; %NaN are used to prevent the patch looping back
    y = [eigF(i,:) NaN ];
    c = squeeze(abs(eigV(:,i,:)))'; c = [ c; NaN NaN NaN ];
    p = patch(x,y,x, 'CDataMapping', 'direct', 'FaceColor', 'interp', 'EdgeColor', 'interp', 'FaceAlpha', 0);
    p.FaceVertexCData = c; p.LineWidth = 2; hold on;
end
hold off;

xlabel('Zeeman splitting (MHz)', 'FontSize', 14);
ylabel('V/h (MHz)', 'FontSize', 14);

set(gcf, 'Color', 'w');
box on;