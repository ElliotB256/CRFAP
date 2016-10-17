%% Test Rabi Freq
% Tests the Rabi freq is correct for F=1 atoms

RF = 2';
Rabi = 0.1';
Bs=(0:0.2:3);

[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RF, Rabi * 2/3, 'iterations', 3, 'F', 2);
plot(B, F, '.');

%%
% Plot Rabi flopping between states to test the MRF.Propagator function.
% Show that after time=Rabi freq we get a revival and revert back to
% initial state.

H = MRF.Hamiltonian(RF, RF, Rabi, 'F', 1);

% Calculate the wavefunction
t = linspace(0, 2*pi/Rabi, 30); t(1) = [];
wavefunc = repmat([ 1 0 0 ]', 1, length(t));
for j=1:length(t);
    tempU = MRF.Propagator(H, t(j));
    wavefunc(:,j) = tempU * wavefunc(:,j);
end

figure(2);
for i=1:3;
    switch i
        case 1
            str = '|mF=+1>';
            fs = [ 1 0 0 ]';
        case 2
            str = '|mF= 0>';
            fs = [ 0 1 0 ]';
        case 3
            str = '|mF=-1>';
            fs = [ 0 0 1 ]';
    end
    fs = repmat(fs, 1, length(t));
    subplot(4, 1, i);
    plot(t*Rabi/2/pi,abs(sum(fs .* wavefunc, 1)).^2);
    title(str);
end

subplot(4, 1, 4);
plot(t*Rabi/2/pi, sum(wavefunc .* conj(wavefunc), 1));
title('Normalisation of wavefunction'); ylim([0 1]);

pos = get(gcf, 'Position');
set(gcf, 'Position', [ 200 200 400 900 ]);

%%
% The same for F=2.
Rabi = 0.5;
H = MRF.Hamiltonian(RF, RF, Rabi * 2/3, 'F', 2);

t = linspace(0, 2*pi/Rabi, 50); t(1) = [];
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

%% 
% Use propagator to examine wavefunction at a time t=1/RF later.
%% 
% First up: F = 1
Rabi = 0.1;
RF = 2;
H = MRF.Hamiltonian(RF, RF, Rabi, 'F', 1);
P = MRF.Propagator(H, 2*pi/Rabi);

% consider action on pure mF states one by one.
a = [ 1 0 0 ]';
pa = conj(a)' * P * a;

b = [ 0 1 0 ]';
pb = conj(b)' * P * b;

c = [ 0 0 1 ]';
pc = conj(c)' * P * c;

fprintf('Probability for pure mF states to return after one Rabi oscillation are:\n')
fprintf('|mF=1 >: %.3f\n', pa)
fprintf('|mF=0 >: %.3f\n', pb)
fprintf('|mF=-1>: %.3f\n', pc)
fprintf('(These numbers should all equal 1.)\n')
assert(abs(pa - 1) < 1e-4, 'Revival expects pa=1');
assert(abs(pb - 1) < 1e-4, 'Revival expects pb=1');
assert(abs(pc - 1) < 1e-4, 'Revival expects pc=1');

%%
% Next up: F = 2
H = MRF.Hamiltonian(RF, RF, Rabi * 2/3, 'F', 2);
P = MRF.Propagator(H, 2*pi/Rabi, 2);

% consider action on pure mF states one by one.
a = [ 1 0 0 0 0 ]';
pa = conj(a)' * P * a;

b = [ 0 1 0 0 0]';
pb = conj(b)' * P * b;

c = [ 0 0 1 0 0 ]';
pc = conj(c)' * P * c;

fprintf('Probability for pure mF states to return after one Rabi oscillation are:\n')
fprintf('|mF=2>: %.3f\n', pa)
fprintf('|mF=1>: %.3f\n', pb)
fprintf('|mF=0>: %.3f\n', pc)
fprintf('(These numbers should all equal 1.)\n')
assert(abs(pa - 1) < 1e-4, 'Revival expects pa=1');
assert(abs(pb - 1) < 1e-4, 'Revival expects pb=1');
assert(abs(pc - 1) < 1e-4, 'Revival expects pc=1');