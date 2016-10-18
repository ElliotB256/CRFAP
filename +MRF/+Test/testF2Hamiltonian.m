%% Test F=2 Hamiltonian
% Tests the F=2 Hamiltonian for specific values with known results.

%%
% Test 1: Separation between bare states.

gFuBB   = 0;
RF     = 1;
Zeeman = 1;
H = MRF.Hamiltonian(Zeeman, RF, gFuBB, 'F', 2);

h = H(0);

% Check we get -Zsf 0 Zsf along diagonal from H0
assert(all(abs(diag(h) - [ 2*Zeeman Zeeman 0 -Zeeman -2*Zeeman ]') < eps), 'H0 should be split by specified Zeeman splitting, with mF basis [ +2 +1 0 -1 -2 ]''');
fprintf('%s: Test 1 passed.\n', mfilename)

%%
% Test 2: Form of matrix for coherence terms.

gFuBB   = 1;
RF     = 1;
Zeeman = 0;

H = MRF.Hamiltonian(Zeeman, RF, gFuBB, 'F', 2, 'theta', 0);
h = H(0);

% get terms we expect to be zero
eZ = [ h(1,1) h(3,1) h(4,1) h(5,1) h(2,2) h(2,4) h(2,5) h(3,1) h(3,3) h(3,5) h(4,1) h(4,2) h(4,4) h(5,1) h(5,2) h(5,3) h(5,5) ];

assert(all(eZ < eps), 'Coherences should be off-diagonal');
fprintf('%s: Test 2 passed.\n', mfilename)

%%
% Test 3: Magnitude of coherence terms
gFuBB  = 1;
RF     = 1;
Zeeman = 0;

H = MRF.Hamiltonian(Zeeman, RF, gFuBB, 'F', 2, 'theta', 0);
h = H(0);

assert(all(abs(abs([ h(1,2) h(2,3) h(3,4) h(4,5) ]) - gFuBB * [ 2 sqrt(6) sqrt(6) 2 ] / 2) < eps), 'Magnitude of coherence terms, Clebsch-Gordan coefficients');
assert(all(abs(abs([ h(2,1) h(3,2) h(4,3) h(5,4) ]) - gFuBB * [ 2 sqrt(6) sqrt(6) 2 ] / 2) < eps), 'Magnitude of coherence terms, Clebsch-Gordan coefficients');
fprintf('%s: Test 3 passed.\n', mfilename)

%%
% Test 4: Periodicity of the Hamiltonian
gFuBB   = 1;
RF     = 1;
Zeeman = 0;
period = 2*pi/RF;

H = MRF.Hamiltonian(Zeeman, RF, gFuBB, 'F', 2, 'theta', 0);
h1 = H(0);
h2 = H(period);

% Render for general intruige
t = linspace(0, period, 100);
Hs = zeros(25, size(1, 1));
for i = 1:length(t);
   h = H(t(i));
   Hs(:, i) = h(:);
end
for i=1:25
    subplot(5,5,i);
    plot(real(Hs(i,:))'); hold on; plot(imag(Hs(i,:))'); hold off;
end

assert(all(h1(:) - h2(:) < eps), 'Hamiltonian must be periodic with periodicity of system (RF)');
fprintf('%s: Test 4 passed.\n', mfilename)

%% 
% Test 5: Hamiltonian is Hermitian

gFuBB   = 1;
RF     = 1;
Zeeman = 0;
period = 2*pi/RF;

H = MRF.Hamiltonian(Zeeman, RF, gFuBB, 'F', 2, 'theta', 0);

ts = linspace(0, period, 20);
for t=ts
   h = H(t);
   assert(all(all(h - conj(h') < eps)), 'Hamiltonian must be a Hermitian matrix');
end
fprintf('%s: Test 5 passed.\n', mfilename)