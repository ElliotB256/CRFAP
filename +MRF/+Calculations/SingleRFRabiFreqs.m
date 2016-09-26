%% Single RF Rabi Freq
% Calculates the single RF rabi frequencies for circular polarised rf in
% the shell trap.

%%
% First load the calculated quadrupole gradient from the 'quadGradient.mat'
% file.
thisloc = fileparts(mfilename('fullpath'));
load([thisloc filesep 'quadGradient.mat']);

%%
% This cell array holds parameters for each dressing rf, including the
% range to calculate the potential over, the rf itself and the resonant
% frequency.

 runs = [ ...
     struct('RF', 3, 'ZeemanSplit', 2.8:0.2:3.8, 'resonance', 3+0.392);
     struct('RF', 3.6, 'ZeemanSplit', 3.0:0.3:4.6, 'resonance', 3.6+0.410);
     struct('RF', 4.2, 'ZeemanSplit', 3.8:0.2:5.0, 'resonance', 4.2+0.390);
     ];

% Data from logbook2016_06_08. Fitted 'by eye' for now.
runs = [ ...
    struct('RF', 3, 'ZeemanSplit', 2.8:0.2:3.8, 'resonance', 3+0.392);
    struct('RF', 3.6, 'ZeemanSplit', 3.0:0.3:4.6, 'resonance', 3.6+0.460);
    struct('RF', 4.2, 'ZeemanSplit', 3.8:0.2:5.0, 'resonance', 4.2+0.410);
    ];

% Data from 29-31 August 2016
runs = [ ...
    struct('RF', 3, 'ZeemanSplit', 2.8:0.2:3.8, 'resonance', 3+0.3760);
    struct('RF', 3.6, 'ZeemanSplit', 3.0:0.3:4.6, 'resonance', 3.6+0.4431);
    struct('RF', 4.2, 'ZeemanSplit', 3.8:0.2:5.0, 'resonance', 4.2+0.4043);
    ];

%%
% Iterate through the runs. For each, calculate the Rabi frequency that
% would result in the measured resonance.
%
% Select the RF transition we are driving with the measured experimental
% resonance. I will plot this in a different color at the end just to
% highlight it.
transitionIndex = 3;
getTrans = @(spectra) spectra(transitionIndex, :);

results = struct('rabi', []);
debug = 1;
Rabis = [0.1:0.05:0.8];


for j=1:length(runs)
    run = runs(j);
    
    if debug
        figure(1);
    end
    
    spectra = [];
    for Rabi=Rabis
        [spec, debugData] = MRF.Spec.Calc(run.ZeemanSplit, run.RF, Rabi, 'qdrpGrad', RFSpecQuad, 'ladderN', 8);
        spectra(:,end+1) = spec;
        
        if debug
            plot(debugData.B, debugData.trapped, '.-', 'Color', [0.5 0.5 0.5]); hold on;
            plot(debugData.B(debugData.min), debugData.trapped(debugData.min), '.', 'Color', [0.8 0.2 0.2]); hold on;
        end
    end
    hold off;
    
    if debug
        xlabel('Bare Zeeman splitting (MHz)');
        ylabel('Dressed state energies (MHz)');
        title('Gravitational sag');
        
        figure(2);
        plot(Rabis, spectra', '-', 'Color', [0.5 0.5 0.5]); hold on
        plot(xlim, run.resonance*[1 1], ':b');
        plot(Rabis, getTrans(spectra), '--', 'Color', [0.8 0.5 0.5], 'LineWidth', 1.5); hold off
        
        xlabel('Rabi Frequency (MHz)');
        ylabel('Transition Energy (MHz)');
        title('RF Spec transitions v Rabi Freq (4 MHz)');
    end
    
    % Now use interp function to find where spec==measured
    st = find(getTrans(spectra) - run.resonance < 0, 1, 'last');
    en = find(getTrans(spectra) - run.resonance > 0, 1, 'first');
    dr = (Rabis(en)-Rabis(st))/1000;
    vals = Rabis(st):dr:Rabis(en);
    
    f = interp1(Rabis, getTrans(spectra), vals, 'linear');
    
    % find the closest agreement over finer range:
    [~,i] = min(abs(f - run.resonance));
    
    if debug
        figure(2);
        hold on; plot([vals(i) vals(i)], ylim, 'b:'); hold off
        pause(0.1);
    end
    
    results(j).rabi = vals(i);
end

%%
% Having calculated the resulting RF Rabi freqs at the '100%' value we
% normalise to, save the results.
thisloc = fileparts(mfilename('fullpath'));
Rabi = [results(:).rabi];
save([thisloc filesep 'rabiFreqs.mat'], 'Rabi');