%% MRF Well Positions
% This calculates the position of wells in the MRF scheme. The quadrupole
% gradient is used to add in the correct offsets with respect to gravity

%% Load Calibrations
% Load calibrations of the quadrupole gradient and magnification
thisLoc = mfilename('fullpath');

quadCalFile = [fileparts(thisLoc) filesep 'quadGradient.mat'];
if ~exist(quadCalFile, 'file')
    error('run QuadGradient_top.m first!');
end
load(quadCalFile);

magFile = [fileparts(thisLoc) filesep 'magnification.mat'];
if ~exist(magFile, 'file')
    error('run Magnification.m first!');
end
load(magFile);

rabiFreqFile = [fileparts(thisLoc) filesep 'rabiFreqs.mat'];
if ~exist(rabiFreqFile, 'file')
    error('run SingleRFRabiFreqs.m first!');
end
load(rabiFreqFile);

clear quadCalFile magFile rabiFreqFile
Rabi = Rabi';

%%
% Define the experimental setup:

Bs = 2.5:0.05:4.5;
Factor = 0.98*[ 0.5 1 1.1 ]';
RFs = [3.0 3.6 4.2 ]';
BarrierHeights = 0.1:0.1:1.6;

% Was having an issue with the network drive and Matlab not finding files.
% Using a local copy of data now.
experimentalResultsFolder = 'C:\Data\MRF\DoubleWells';


%%
% Enumerate through all barrier heights

minima = cell(length(BarrierHeights), 1); % cell required as number of minima may change

for i=1:length(BarrierHeights)
    fprintf('Generating MRF potential (%d of %d) ...\n', i, length(BarrierHeights))
    
    Barrier = BarrierHeights(i);
    Factor(2) = Barrier * 0.98;
    ExptRabi = Rabi .* Factor;
    
    [ Fs, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, ExptRabi, 'iterations', 4, 'qdrpGrad', QuadGradDblShell);
    Fs = MRF.sortEnergies(B, MRF.ladder(RFs, 20, Fs));  
    
    % include gravitational sag
    Fs = Fs + repmat(MRF.gpe(B, QuadGradDblShell), size(Fs,1), 1);
        
    % Extract minima from the MRF Levels
    minima{i} = MRF.Calculations.getMinima(B, Fs(3,:));
    
end

fprintf('Loop complete\n')

%%
% Sort resulting minima into two position lists

minPos = zeros(length(minima),2);
for i=1:length(minima)
    if length(minima{i}) == 1
        minPos(i, 2) = minima{i};
        minPos(i, 1) = NaN;
    elseif isempty(minima{i})
        minPos(i, :) = [NaN NaN];
    else
        minPos(i, :) = minima{i};
    end
end

%%
% Now we must load the experimental data for comparison purposes
currD = cd;
cd(experimentalResultsFolder);
shell_densities_onlyBoth

clear shots loadImg absImgF roi rectx recty;
clear mag pixelSize quadC;
clear pix2um;

cla;
striplinePlot(bhs, navgcys, fliplr(exptY))

cd(currD);


%%
% Convert these minima into actual positions using the quadrupole gradient.

minPosum = minPos ./ QuadGradDblShell * 1e4 -1000;

hold on; plot(BarrierHeights, minPosum, 'LineWidth', 3, 'Color', [0.2 0.2 0.2]); hold off;



