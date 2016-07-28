%% MRF Well Positions
% This calculates the position of wells in the MRF scheme. The quadrupole
% gradient is used to add in the correct offsets with respect to gravity

%%
% Define the experimental setup:

Bs = 2.5:0.1:4.5;
BaseRabi = [ 338.7 401.2 357.0]' * 1e-3;
RFs = [3.0 3.6 4.2 ]';
qdrpGrad = 62.4511*0.96; % Gauss/cm at 20 A
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
    Rabi = BaseRabi; Rabi(2) = Barrier * BaseRabi(2);
    
    [ Fs, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', 4, 'qdrpGrad', qdrpGrad);
    Fs = MRF.sortEnergies(B, MRF.ladder(RFs, 20, Fs));  
        
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

minPosum = -minPos ./ qdrpGrad * 1e4 + 250;

hold on; plot(BarrierHeights, minPosum, 'LineWidth', 3, 'Color', [0.2 0.2 0.2]); hold off;



