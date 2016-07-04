%% Testing out APs in matlab

freq = [];

% iterate over all magnetic fields
Bs=0:0.025:2;
for B=Bs;
    
    % Define the Hamiltonian. RF, Rabi
    RF = 1; % MHz
    Rabi = 0.2; % Rabi freq in MHz
    H = MRF.Hamiltonian(B, RF, Rabi);
    period = 1/RF;
    
    % time integration
    dt = 0.01 * period;
    DM = eye(3);
    for t=0:dt:period
        h = H(t);
        deltaM = -1i * dt * ( h * DM );
        DM = DM + deltaM;
    end
    e = eig(DM);
    freq(:,end+1) = angle(e)/period;
end

freq = sort(freq, 1);
plot(Bs, freq, 'k.');%; freq+2*period; freq+3*period; freq+4*period], '.k')

%%
% Check propagator against that from Mathematica
RFs = [3 3.6 4.2 ]';
BRFs = 0.8 * [ 0.5 0.2 1.1 ]';
H = MRF.Hamiltonian(1, RFs, BRFs);
U = MRF.propagator(H, 10*pi/3);
e = eig(U);
eigF = angle(e)/(10*pi/3)