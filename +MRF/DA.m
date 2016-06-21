%% Pissing about with dressed atoms

B = 3;
RF = 6;
rabi = pi/4;
H = @(t) [ -B 0 0; 0 B 0; 0 0 B ] + [ 0 rabi*cos(RF*t) 0; rabi*cos(RF*t) 0 rabi*cos(RF*t) ; 0 rabi*cos(RF*t) 0 ];

psi0 = [ 1 0 0 ]';

% integrate psi as a function of time
cpsi = psi0;
dt = 0.0001;
t=0:dt:5;
psi = zeros(3, length(t));
for j=1:length(t)
    ct = t(j);
    cpsi = cpsi ./ (sum(conj(cpsi) .* cpsi, 1))^0.5;
    psi(:,j) = cpsi;
    
    cpsi = cpsi - 1i * dt * H(ct) * cpsi;
end
plot(t, conj(psi).*psi);

%%
% Plot difference between mF states as a function of time
a = psi(1,:);
b = psi(2,:);
plot(t, abs(conj(a).*a - conj(b).*b));