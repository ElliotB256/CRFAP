% function [ allMinimaCM ] = mrfCOMs(amp30, amp36, amp42)
% Get spatial positions of trap COMs in um from quad centre given mrf 
% amplitudes 
% (coefficient * 100% barrier heights)
% Only accounts for ramping barrier - needs constant f1 and f3 - might need
% to change this in future!

%% Trap centre spatial positions

%%
% First define trap parameters

RFs = [3.0 3.6 4.2 ]';

rabi100 = [0.391 0.459 0.410]'; % '100%' reference values measured in single-rf traps (MHz)

c1 = 0.5*rabi100(1);
c2 = [0.20]*rabi100(2); % the final barrier heights
c3 = 1*rabi100(3);

% Work with series 2 for now
amp30 = c1 * rabi100(1);
amp36 = c2 * rabi100(2); % the final barrier heights
amp42 = c3 * rabi100(3);

qdrpGrad = 62.4511; % Gauss/cm at 20 A 
Bs=(2.5:0.01:5); % range of Bs to consider 


%% Trap minima

% trapMinima = locateWells(Bs, FworkingI, 8)

% TwoMinList = [];
% OneMinList = [];
ladno = 3;

count2 = 0;
count1 = 0;
% TwoMinList = [];
% OneMinList = [];
for in = 1:length(amp36)
    rabis = [amp30 amp36(in) amp42]';
    disp(amp36(in)/0.357)
    
    [ Fs, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, rabis, 'iterations', 4, 'qdrpGrad', qdrpGrad);
    Fs = MRF.sortEnergies(B, MRF.ladder(RFs, ladno, Fs));  
    FsTest = Fs(1,:);
    
    Fws = Fs(1,:)+repmat(MRF.gpe(B, qdrpGrad), size(Fs(1,:),1), 1);
    FwI(:,:,in) = interp1(B, Fws, Bs);
    
%     amp36(in)
%     locateWells(Bs, FwI(:,:,in), 8)
    trapMinTemp = MRF.locateWells(Bs, FwI(:,:,in), 8);
    if size(trapMinTemp,1) == 2
        count2 = count2+1;
%         TwoMinList = [TwoMinList locateWells(Bs, FwI(:,:,in), 8)]
        TwoMinList(:,2:3,count2) = MRF.locateWells(Bs, FwI(:,:,in), 8);
        TwoMinList(:,1,count2) = amp36(in);
    elseif size(trapMinTemp,1) == 1
        count1 = count1+1;
%         OneMinList = [OneMinList locateWells(Bs, FwI(:,:,in), 8)]
        OneMinList(:,2:3,count1) = MRF.locateWells(Bs, FwI(:,:,in), 8);
        OneMinList(:,1,count1) = amp36(in);
    else 
        disp('How many minima do we have??')
    end
  
end

%%
% rearrange data 
count2;
count1;
if count2 == 0
    allMinima = reshape(permute(OneMinList, [1 3 2]), size(OneMinList,1)*size(OneMinList,3), size(OneMinList,2));
elseif count1  == 0
    reshape(permute(TwoMinList, [1 3 2]), size(TwoMinList,1)*size(TwoMinList,3), size(TwoMinList,2));
else 
    allMinima = [reshape(permute(TwoMinList, [1 3 2]), size(TwoMinList,1)*size(TwoMinList,3), size(TwoMinList,2)); ...
    reshape(permute(OneMinList, [1 3 2]), size(OneMinList,1)*size(OneMinList,3), size(OneMinList,2))];
end


% convert B in G to separation in um
% have 3.1226 G/cm/A; working at 20A so 62.452 G/cm
allMinimaCM = [allMinima allMinima(:,2,:)/(qdrpGrad/1e4)];




% end
