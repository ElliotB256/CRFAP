%% Small Frequency separation MRF Spectroscopy
RF = [ 4.0 4.1 4.2 ]';
%Bs =  3.9:0.025:4.25; %4.05:0.02:4.25; %3.8:0.025:4.4;
Bs = 3.9:0.05:4.3;
Rabi = [ 134 144 1.1*136 ]' * 1e-3;
BaseRabi = Rabi;
BarrierPct  = (0.4 : 0.025: 1.2)*0.9;
QdrpGrad = 10*62.5*0.96;

spectra = [];

%bloody trapped detection not working nicely
BsvBarrier = { ...
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.15 4.16 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    [ 4.0 4.01 4.24 4.25 ];
    };
    
    

for i=1:length(BarrierPct)%BarrierPct
    bpct = BarrierPct(i);
    Bs = BsvBarrier{i};
    fprintf('Barrier at %.0f percent.\n', bpct*100)
    Rabi(2) = bpct * BaseRabi(2);
    [spec, debug] = MRF.Spec.Calc(Bs, RF, Rabi, 'ladderN', 200, 'qdrpGrad', QdrpGrad, ...
        'minMethod', 'settle', 'minStart', 4.15);
    spectra(:,end+1) = spec;
    
    plot(debug.B, debug.trapped, '.-', 'Color', [0.5 0.5 0.5]); hold on;
    plot(debug.B(debug.min), debug.trapped(debug.min), '.', 'Color', [0.8 0.2 0.2], 'MarkerSize', 10); hold on;
end
hold off;
title('MRF APs used for rf spectroscopy calc');
xlabel('Bare Zeeman Splitting (MHz)');
ylabel('Energy (MHz)');

figure(2);
plot(BarrierPct * 100, spectra', '.-', 'Color', [0.8 0.2 0.2]); ylim(4.2+[0.1 0.4]);

xlabel('Barrier Height (% of flat potential value)');
ylabel('Transition Energy (MHz)');
title('RF Spec transitions v Rabi Freq');