%% MRF Shifts
% Calculate the shift in resonance position as a function of the RF
% amplitudes. Compare this to first order theory estimate.

RFs = [ 3 3.6 ]';
gFuBB = [ 0.1 0.1 ]';
qdrpGrad = 100;

Zsf=(2.5:0.2:3.6);
[ F, B ] = MRF.MeshedQuasiEnergies(Zsf, RFs, gFuBB, 'iterations', 7, 'qdrpGrad', qdrpGrad, 'F', 1);
F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 30, F));

plot(B,F2)


%%
% Step through different amplitudes of \Omega_2 and observe shift in
% resonance position of resonance near \omega_1.

rfAmps = 0.05:0.05:0.4;

minBs = zeros(1, length(rfAmps));

for i=1:length(rfAmps)
    a = rfAmps(i);
    
    [ F, B ] = MRF.MeshedQuasiEnergies(Zsf, RFs, [gFuBB(1) a]', 'iterations', 7, 'qdrpGrad', qdrpGrad, 'F', 1);
    F2 = MRF.sortEnergies(B, MRF.ladder(RFs, 30, F));
    
    % get minima of energy levels
    [~,j] = min(F2(3,:));
    minBs(i) = B(j);
    
    plot(B, F2(3,:), '-'); hold on; plot(B(j), F2(3,j), '.', 'MarkerSize', 12); hold off; pause(0.01);
end

%%
% Plot the shift in the resonance arising from the second rf.

shift = (0.16 / 0.26)* (rfAmps .^ 2) / (RFs(2)-RFs(1));

plot( rfAmps, minBs - 3.0);
hold on; plot( rfAmps, shift); hold off;
