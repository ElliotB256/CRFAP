%%
% Test the periodicity of the Hamiltonian is correct
B = 3;
RF = 3;
Rabi = 1;

H = MRF.Hamiltonian(B, RF, Rabi, 'F', 1);

period = 2*pi/RF;
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

% Should see the Hamiltonian periodic in 1/RF.