
% Express the wavefunction in the mF basis. F=1 system.


BRFs = 1; RFs = 3;

% todo: set period automatically from RFs
%fund = MRF.GetFundamental(RFs);
fund = 3;
freq = [];

Bs=0:0.025:5;
%Bs = Bs/1000;
for B=Bs;
H = MRF.Hamiltonian(B, RFs, BRFs);

% time integration
dt= 0.01/fund;
DM = eye(3);
dbg = [];
for t=0:dt:(1/fund)
    h = H(t);
    deltaM = -1i * dt * (h * DM );
    DM = DM + deltaM;
    dbg(end+1) = [ 1 0 0 ] * DM * [ 1 0 0 ]';
end
e = eig(DM);
freq(:,end+1) = angle(e)*fund;
end

freq = sort(freq, 1);
plot(Bs, freq, 'k.');%; freq+2*period; freq+3*period; freq+4*period], '.k')


%%
plot(conj(dbg) .* dbg)

%% Plot hamiltonian as a function of time
H = MRF.Hamiltonian(3*0.7, RFs,5* BRFs);
ts=0:dt:(1/fund);
hs = zeros(9,length(ts));
for i=1:length(ts)
    t = ts(i);
    h = H(t);
    hs(:,i) = h(:);
end

plot(ts,hs')