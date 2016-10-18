%% Gravitational Sag
% This graph attempts to demonstrate how two separate species can be
% independently controlled and manipulated in a 2-rf dressed trap provided
% that their gF factors differ. We take the example of Rb85 and Rb87, whose
% gF factors differ by 2/3.
%
% This graph aims to show how by controlling one of the Rabi freqs we can
% raise and lower one of the species through the well formed of the other
% species.

%% Section 1: Calculate trap paramters for Rb87 single condensate well.
% I want to calculate the Thomas-Fermi radius of an 87 condensate in a well
% so that we can show we can move the impurity from one side of it to the
% other just using the Rabi frequency. I calculate the trap frequencies
% using the single rf case, which is ok for the illustrative purposes
% described here.

RF1 = 4.5;
QdrpGrad = 100;
RF1A = 0.4;

w = SRF.shellTrapFrequencies(RF1, RF1A / 0.7, QdrpGrad);

a       = 100 * Constants.bohr;
N       = 1e5;            % 10 thousand atoms
omega   = geomean(w);          % geometric trap frequency, Hz

oscLength = @(m, angFreq) (Constants.hbar ./ (Constants.amu * m * angFreq)).^0.5;

chemPot = 1/2* ( 15 * N * a / oscLength(87, omega)) ^ (2/5) * omega;

% Pethick p155
TFradius = 1e6 * (2 * (chemPot * Constants.hbar) ./ (87 * Constants.amu) ./ ((2 * pi * w(3)).^2)).^0.5;

fprintf('87 Condensate parameters:\n')
fprintf('Chemical potential (kHz): %.2f\n', chemPot/(2 * pi * 1000))
fprintf('Atom number\t\t\t\t: %.1e\n', N)
fprintf('Trap frequencies (Hz)\t: (%.2f, %.2f, %.2f)\n', w)
fprintf('T-F radius (um)\t\t\t: %.2f\n', TFradius)

%% Section 2: Calculate trap frequencies for 85
% Calculate the trap frequencies for Rb85 over a range of omega. Use the
% two frequency calculation here. Determine the gravitational sag in the
% potential. One trap frequencies are found we can use this to determine
% the values of oscillator lengths

RF2  = 3;
amps = 0.05:0.05:0.4;
amps = 0.2:0.02:0.3;
amps = 0.1:0.02:0.2;
uB   = 2.9:0.2:3.4;
RFs = [RF1 RF2]';

fz = ones(1, length(amps));
zPos85 = ones(1, length(amps));
zPos85_unshifted = ones(1, length(amps));

for i=1:length(amps)
    RF2A = amps(i);
    
    gF = 0.7 * 2/3;
    mass = 85;
    
    % NOTE: the factor of 2/3 on the rf amplitudes is the ratio of gF
    % factors.
    [ F, zsfs ] = MRF.MeshedQuasiEnergies(uB, RFs, 2/3*[RF1A RF2A]', 'iterations', 9, 'qdrpGrad', QdrpGrad, 'gF', gF, 'mass', mass, 'F', 2);
    F2 = MRF.sortEnergies(zsfs, MRF.ladder(RFs, 10, F), 'F', 2);
    lad = MRF.ladder(RFs, 3, F2);
    glad = lad + repmat(MRF.gpe(zsfs, QdrpGrad, 'gF', gF, 'mass', mass), size(lad,1), 1);
    
    % select the trapped manifold
    trapped = glad(2, :);
    
    % calculate the trap frequency (Note: in Hz!)
    z = zsfs / gF / QdrpGrad * 1e4;
    
    % calculate the minimum position
    [~,j] = min(trapped);
    zPos85(i) = z(j);
    
    % Note: only take a small section around the minimum for trap freq fit
    mask = abs(z-z(j)) < 8;
    
    fs = Util.getTrapFreq(z(mask), trapped(mask), mass);
    fprintf('Trap freq fit = %.2f Hz, uncertainty = %.1f%%\n', fs(2), (fs(3)-fs(2))/fs(2)*100)
    fz(i) = fs(2);
   
    
    plot(zsfs, glad, '.-');
    hold on; plot(zsfs(j), trapped(j), '.', 'MarkerSize', 10); hold off;
    pause(0.1);
    
    % For comparative purposes - find the trap minimum in absence of the
    % shift caused by the other dressing RF component.
    [ F, zsfs ] = MRF.MeshedQuasiEnergies(uB, RF2, RF2A', 'iterations', 8, 'qdrpGrad', QdrpGrad, 'gF', gF, 'mass', mass, 'F', 2);
    F2 = MRF.sortEnergies(zsfs, MRF.ladder(RFs, 10, F), 'F', 2);
    lad = MRF.ladder(RFs, 3, F2);
    glad = lad + repmat(MRF.gpe(zsfs, QdrpGrad, 'gF', gF, 'mass', mass), size(lad,1), 1);
    trapped = glad(2,:);
    [~,j] = min(trapped);
    z = zsfs / gF / QdrpGrad * 1e4;
    zPos85_unshifted(i) = z(j);
    
end

%% Section 3: Calculate harmonic oscillator lengths for 85 at each Rabi freq

a85 = oscLength(85, 2*pi*fz) * 1e6;

%% Section 4: Determine gravitational sag for 87 in MRF potential

uB = 4.0:0.1:4.8;
zPos87 = ones(1, length(amps));
zPos87_unshifted = ones(1, length(amps));

for i=1:length(amps)
    RF2A = amps(i);
    
    gF = 0.7;
    mass = 87;
    
    [ F, zsfs ] = MRF.MeshedQuasiEnergies(uB, RFs, [RF1A RF2A]', 'iterations', 9, 'qdrpGrad', QdrpGrad, 'gF', gF, 'mass', mass);
    F2 = MRF.sortEnergies(zsfs, MRF.ladder(RFs, 10, F));
    lad = MRF.ladder(RFs, 3, F2);
    glad = lad + repmat(MRF.gpe(zsfs, QdrpGrad, 'gF', gF, 'mass', mass), size(lad,1), 1);
    
    % select the trapped manifold
    trapped = glad(2, :);
    
    % calculate the trap frequency (Note: in Hz!)
    z = zsfs / gF / QdrpGrad * 1e4;
    
    % calculate the minimum position
    [~,j] = min(trapped);
    zPos87(i) = z(j);
    
    plot(zsfs, glad, '.-');
    hold on; plot(zsfs(j), trapped(j), '.', 'MarkerSize', 10); hold off;
    pause(0.1);
    
    % For comparative purposes - find the trap minimum in absence of the
    % shift caused by the other dressing RF component.
    [ F, zsfs ] = MRF.MeshedQuasiEnergies(uB, RF1, RF1A', 'iterations', 6, 'qdrpGrad', QdrpGrad, 'gF', gF, 'mass', mass);
    F2 = MRF.sortEnergies(zsfs, MRF.ladder(RF1', 10, F));
    lad = MRF.ladder(RF1', 3, F2);
    glad = lad + repmat(MRF.gpe(zsfs, QdrpGrad, 'gF', gF, 'mass', mass), size(lad,1), 1);
    trapped = glad(2,:);
    [~,j] = min(trapped);
    z = zsfs / gF / QdrpGrad * 1e4;
    zPos87_unshifted(i) = z(j);
    
end

%% Assemble the graph.
expectedPos = RFs(1) / 0.7 / QdrpGrad * 1e4; 
plot(amps,  zPos87-expectedPos, 'r-'); hold on;
plot(amps,  zPos85-expectedPos, 'b-'); hold on;
plot(amps,  zPos87_unshifted-expectedPos, 'r:'); hold on;
plot(amps,  zPos85_unshifted-expectedPos, 'b:'); hold off;
legend('87', '85');

title('Gravitational sag in the 2-rf trap');
xlabel('3.0 MHz amplitude, $\Omega_1$ (MHz)', 'Interpreter', 'Latex');
ylabel('Gravitational sag ($\mu$m)', 'Interpreter', 'Latex');

%% Second graph version
% This version incorporates the harmonic oscillator lengths and TF radius
% previously calculated.
expectedPos = RFs(1) / 0.7 / QdrpGrad * 1e4; 
figure(2);cla;

% define colors:
c85 = [ 0.8 0.2 0.2 ];
c87 = [ 0.2 0.2 0.8 ];

% draw filled polygon showing the 85 +- half oscillator length.
sag85 = zPos85-expectedPos;
fill([amps fliplr(amps)]*1e3, [sag85+a85/2 fliplr(sag85-a85/2)], c85, 'LineStyle', 'none', 'FaceAlpha', 0.5);
hold on;

% draw a filled polygon representing the 87 position +- TF radius/2
sag87 = zPos87-expectedPos;
fill([amps fliplr(amps)]*1e3, [sag87+TFradius/2 fliplr(sag87-TFradius/2)], c87, 'LineStyle', 'none', 'FaceAlpha', 0.5);
hold off;

% Lines depicting the location of potential minima
hold on;
h1 = plot(amps*1e3, sag87, 'Color', c87);
h2 = plot(amps*1e3, sag85, 'Color', c85);
plot(amps*1e3, zPos87_unshifted-expectedPos, '--', 'Color', c87);
plot(amps*1e3, zPos85_unshifted-expectedPos, '--', 'Color', c85);
hold off

%xlim([250 400]);
xlabel('4.5 MHz RF amplitude $\Omega_1$ (kHz)', 'Interpreter', 'Latex');
ylabel('Gravitational sag ($\mu$m)', 'Interpreter', 'Latex');

legend([h1 h2], {'87', '85'}, 'Location', 'SouthEast')

%% Third Graph: Demonstrate the APs for both species on same axis
% 

B = 1:0.25:7;

% Calculate the potential for 87:
gF = 0.7;
mass = 87;
[ F, zsfs ] = MRF.MeshedQuasiEnergies(B, RFs, [RF1A RF2A]', 'iterations', 7, 'qdrpGrad', QdrpGrad, 'gF', gF, 'mass', mass);
E87 = MRF.sortEnergies(zsfs, MRF.ladder(RFs, 10, F));
lad87 = MRF.ladder(RFs, 3, E87);
glad87 = lad87 + repmat(MRF.gpe(zsfs, QdrpGrad, 'gF', gF, 'mass', mass), size(lad87,1), 1);
z87 = zsfs / (QdrpGrad * gF);


% Same for 85
gF = 0.7 * 2/3;
mass = 85;
[ F, zsfs ] = MRF.MeshedQuasiEnergies(B, RFs, [RF1A RF2A]', 'iterations', 7, 'qdrpGrad', QdrpGrad, 'gF', gF, 'mass', mass, 'F', 2);
E85 = MRF.sortEnergies(zsfs, MRF.ladder(RFs, 10, F), 'F', 2);
lad85 = MRF.ladder(RFs, 3, E85);
glad85 = lad85 + repmat(MRF.gpe(zsfs, QdrpGrad, 'gF', gF, 'mass', mass), size(lad85,1), 1);
z85 = zsfs / (QdrpGrad * gF);

plot(z85, glad85, 'Color', c85);
hold on;
plot(z87, glad87, 'Color', c87);
hold off;

% This figure is not particularly useful...