%%
% Test the AP calculation for F=1 atoms

%%
% Test 1: The AP on resonance has a splitting given by the Rabi freq
RFs = 2;
gFuBB = 0.5;
F = MRF.GetQuasiEnergies(RFs, RFs, gFuBB, 'F', 1, 'theta', 0);

levelSep = diff(F);
assert(all(abs(levelSep - gFuBB)/gFuBB < 1e-3), 'Separation of levels on resonance should be Rabi freq for F=1'); 
fprintf('%s: Test 1 passed.\n', mfilename)

%%
% Other possible tests: take a small region around the resonance, check the
% shape is always continous and smooth?