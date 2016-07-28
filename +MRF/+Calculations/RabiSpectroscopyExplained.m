%% Rabi Spectroscopy Explained
% In RF spectroscopy we typically apply an rf probe and attempt to measure
% the location in frequency space of transitions. From this we may
% determine the Rabi frequency of the dressing RF.
%
% To turn it around, in 'Rabi spectroscopy' we vary the Rabi frequency to
% measure the frequency of an unknown probe RF, eg a spurious signal that
% may be limiting the lifetime of the atoms.
%
% Let's take a look at two real sets of data. To take this data, a shell
% was loaded with thermal atoms at ~uK and then held for 8s. The atom
% number was recorded after the hold and the Rabi frequency varied from one
% sequence to the next. Doing this sweeps the base of the trap into and out
% of resonance with the applied spurious rf.
%
% Two such data, for 4.1 MHz and 4.2 MHz, are shown below. The Rabi
% frequencies were calibrated using rf spectroscopy for a single value of
% Rabi frequency.

figure(1);
[~,i] = sort(data41(:,1));
plot(data41(i,1)*BaseRabi41, data41(i,2),'.-');
hold on;
[~,i] = sort(data40(:,1));
plot(data40(i,1)*BaseRabi40, 1e6*data40(i,2), '.-');
hold off;
legend('4.1 MHz', '4.0 MHz');
xlabel('Dressing RF Rabi frequency (kHz)');
ylabel('Remaining atom number');


%%
% The weak probe rf that drive transitions from trapped states causes a
% maximum loss rate (and minimum atom number) when resonant with
% transitions at n\omega+-m\Omega, where \Omega = rabi freq, m = { 1, 2 }
% and \omega = dressing rf frequency. There are some shifts due to gravity
% but we can include these in the numerics.
%
% Through numerics we can calculate the energy splittings (and thus
% expected resonances) for each RF frequency over a range of rabi
% frequencies. We then look for coincidence in the resonant frequencies for
% each data set to determine at what frequency the external probe rf is at.
% It's a one-to-many mapping for each dressing Rabi frequency/atom number
% pair, because that might correspond to a resonance at any value of m=1,2
% or n=0,1,2,3.... This causes a repeating motif for each dressing rf
% frequency. This motif repeats at a frequency equal to the dressing rf.
% Thus, for different dressing rfs the motifs dephase with respect to each
% other. Dips in the signal will be coincident for all plotted dressing rfs
% only at the right value of n. This will help us narrow down our hunt for
% noise.
% 
% Let's plot the resulting graphs:

figure(1);
plot(freqs40, 1e6*atomN40); hold on;
plot(freqs41, atomN41); hold off;
xlabel('Expected resonant frequencies of trap (MHz)');
ylabel('Atom number remaining');
xlim([0 5]);

%% 
% The pattern is badly dephased at n=1 but in phase at n=0. This suggests
% noise on the apparatus may be present at either 350, 550 or 700 kHz
% (though we cannot from the data determine exactly which of these it is).
% The dip at 550 might also be just due to LZ losses for the lower rabi
% freqs or the radial trap frequency vanishing for the higher freqs. Here's
% a closeup of the interesting region:
figure(1);
plot(freqs40, 1e6*atomN40); hold on;
plot(freqs41, atomN41); hold off;
xlabel('Expected resonant frequencies of trap (MHz)');
ylabel('Atom number remaining');
xlim([0.2 0.85]);