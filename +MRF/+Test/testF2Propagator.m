%% Test F=2 propagator
% Runs a series of tests to verify the integrity of the F=2 propagator.

%%
% Test 1: Examine accuracy of collapse and revival mechanics.

gFuBB = 0.1;
RF = 2;
periodicTime = 2*pi/gFuBB;

H = MRF.Hamiltonian(RF, RF, gFuBB, 'F', 2);
P = MRF.Propagator(H, 2*pi/gFuBB, 2);

% consider action on pure mF states one by one.
a = [ 1 0 0 0 0 ]';
pa = conj(a)' * P * a;

b = [ 0 1 0 0 0 ]';
pb = conj(b)' * P * b;

c = [ 0 0 1 0 0 ]';
pc = conj(c)' * P * c;

fprintf('Probability for pure mF states to return after one Rabi oscillation are:\n')
fprintf('|mF=2>: %.3f\n', pa)
fprintf('|mF=1>: %.3f\n', pb)
fprintf('|mF=0>: %.3f\n', pc)
fprintf('(These numbers should all equal 1.)\n')
fprintf('Note: slight discrepancy can be improved by increasing accuracy of MRF.Propagator evaluation.\n');
assert(abs(pa - 1) < 1e-3, 'Revival expects pa=1');
assert(abs(pb - 1) < 1e-3, 'Revival expects pb=1');
assert(abs(pc - 1) < 1e-3, 'Revival expects pc=1');
fprintf('%s: Test 1 passed.\n', mfilename)

%%
% Test 2: Revival does not occur at half this duration
P = MRF.Propagator(H, periodicTime/2, 2);

pa = conj(a)' * P * a;
pb = conj(b)' * P * b;
assert(abs(abs(pa) - 1) > 1e-2, 'Revival should not have occured');
assert(abs(abs(pb) - 1) > 1e-2, 'Revival should not have occured');
fprintf('%s: Test 2 passed.\n', mfilename)

%%
% General intruige: Plot the collapse and revival
% Calculate the wavefunction
Rabi = 0.5;
H = MRF.Hamiltonian(RF, RF, Rabi * 2/3, 'F', 2);

t = linspace(0, 2*pi/(Rabi * 2/3), 50); t(1) = [];
wavefunc = repmat([ 1 0 0 0 0 ]', 1, length(t));
for j=1:length(t);
    tempU = MRF.Propagator(H, t(j), 2);
    wavefunc(:,j) = tempU * wavefunc(:,j);
end

figure(2);
for i=1:5;
    switch i
        case 1
            str = '|mF=+2>';
            fs = [ 1 0 0 0 0 ]';
        case 2
            str = '|mF=+1>';
            fs = [ 0 1 0 0 0 ]';
        case 3
            str = '|mF=-1>';
            fs = [ 0 0 1 0 0 ]';
        case 4
            str = '|mF=-1>';
            fs = [ 0 0 0 1 0 ]';
        case 5
            str = '|mF=-2>';
            fs = [ 0 0 0 0 1 ]';
    end
    fs = repmat(fs, 1, length(t));
    subplot(6, 1, i);
    plot(t*Rabi/2/pi,abs(sum(fs .* wavefunc, 1)).^2);
    title(str);
    ylim([0 1]);
end

subplot(6, 1, 6);
plot(t*Rabi/2/pi, sum(wavefunc .* conj(wavefunc), 1));
title('Normalisation of wavefunction'); ylim([0 1]);

pos = get(gcf, 'Position');
set(gcf, 'Position', [ 200 200 400 900 ]);