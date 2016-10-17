%% Test F=1 Hamiltonian
% Tests the F=1 Hamiltonian for specific values with known results.

%%
% Test 1: Separation between bare states.

Rabi   = 0;
RF     = 1;
Zeeman = 1;
H = MRF.Hamiltonian(Zeeman, RF, Rabi, 'F', 1);

h = H(0);

% Check we get -Zsf 0 Zsf along diagonal from H0
assert(all(abs(diag(h) - [ Zeeman 0 -Zeeman ]') < eps), 'H0 should be split by specified Zeeman splitting');
fprintf('%s: Test 1 passed.\n', mfilename)

%%
% Test 2: Form of matrix for coherence terms.

Rabi   = 1;
RF     = 1;
Zeeman = 0;

H = MRF.Hamiltonian(Zeeman, RF, Rabi, 'F', 1, 'theta', 0);
h = H(0);

assert(all([ h(1) h(3) h(5) h(7) h(9) ] < eps), 'Coherences should be off-diagonal');
fprintf('%s: Test 2 passed.\n', mfilename)

%%
% Test 3: Magnitude of coherence terms
Rabi   = 1;
RF     = 1;
Zeeman = 0;

H = MRF.Hamiltonian(Zeeman, RF, Rabi, 'F', 1, 'theta', 0);
h = H(0);

assert(all(abs([ h(2) h(4) h(6) h(8) ] - Rabi / sqrt(2)) < eps), 'Magnitude of coherence should be gF uB Bi / sqrt(2)');
fprintf('%s: Test 3 passed.\n', mfilename)

%%
% Test 4: Periodicity of the Hamiltonian
Rabi   = 1;
RF     = 1;
Zeeman = 0;
period = 2*pi/RF;

H = MRF.Hamiltonian(Zeeman, RF, Rabi, 'F', 1, 'theta', 0);
h1 = H(0);
h2 = H(period);

% Render for general intruige
t = linspace(0, period, 100);
Hs = zeros(9, size(1, 1));
for i = 1:length(t);
   h = H(t(i));
   Hs(:, i) = h(:);
end
for i=1:9
    subplot(3,3,i);
    plot(real(Hs(i,:))'); hold on; plot(imag(Hs(i,:))'); hold off;
end

assert(all(h1(:) - h2(:) < eps), 'Hamiltonian must be periodic with periodicity of system (RF)');
fprintf('%s: Test 4 passed.\n', mfilename)