function [ allMinimaCM ] = mrfPositions(amp30, amp36, amp42)
% Get trap centre spatial positions in um from quad centre given mrf 
% amplitudes 
% (coefficient * 100% barrier heights)
% Only accounts for ramping barrier - needs constant f1 and f3 - might need
% to change this in future!

%% Trap centre spatial positions

%%
% First define trap parameters

RFs = [3.0 3.6 4.2 ]';

% rabi100 = [0.391 0.459 0.410]'; % '100%' reference values measured in single-rf traps (MHz)

% Work with series 2 for now
% amp30 = c1 * rabi100(1);
% amp36 = c2 * rabi100(2); % the final barrier heights
% amp42 = c3 * rabi100(3);

qdrpGrad = 62.4511; % Gauss/cm at 20 A 
Bs=(2.5:0.01:5); % range of Bs to consider 

%%
% Then calculate potentials. Calculate energies, make a ladder, sort
% energies, then re-ladder to get correct states. Add in gpe. 
% Just one set of parameters!

% % % ladno = 3;
% % % 
% % % rabis = [amp30 amp36(13) amp42]';
% % % [ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, rabis, 'iterations', 4, 'qdrpGrad', qdrpGrad);
% % % F = MRF.sortEnergies(B, MRF.ladder(RFs, ladno, F));
% % % 
% % % figure(5)
% % % plot(B,F)
% % % 
% % % lad = MRF.ladder(RFs, ladno, F);
% % % 
% % % figure(6)
% % % plot(B, lad,'.')
% % % 
% % % figure(1)
% % % plot(B, lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1), '.');
% % % 
% % % Ftest = F(1,:);
% % % 
% % % Fworking = Ftest+repmat(MRF.gpe(B, qdrpGrad), size(Ftest,1), 1);
% % % 
% % % FworkingI = interp1(B, Fworking, Bs); % so this gives me equal spacing to plot lots of them!
% % % 
% % % figure(2)
% % % plot(Bs, FworkingI,'.')
% % %  
% % % % testfunction = interp1(B, Fworking, 'linear', 'pp');
% % % % [breaks,coefs,l,k,d] = unmkpp(testfunction)
% % % % test2 = fminsearch(testfunction, 3.5)
% % % 
% % % % test2 = fminsearch(Fworking, 3.5)

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



%% plot!
% figure(106)
% plot(allMinimaCM(:,1)/rabi100(2), allMinimaCM(:,4),'.')


%%
% look for maxima to see where atoms would leak out - and use this as the
% cutoff point (?)

%%
% for interest, make a corresponding stripline plot of the theory 

% % % % % tempBs = 3:0.02:4.7;
% % % % % 
% % % % % for i = 1:length(amp36)
% % % % %     rabis = [amp30 amp36(i) amp42]';
% % % % %     
% % % % %     [ Fs, B ] = MRF.MeshedQuasiEnergies(tempBs, RFs, rabis, 'iterations', 4, 'qdrpGrad', qdrpGrad);
% % % % %     Fs = MRF.sortEnergies(B, MRF.ladder(RFs, ladno, Fs));
% % % % %     
% % % % % %     lads = MRF.ladder(RFs, ladno, Fs);
% % % % %     
% % % % %     FsTest = Fs(1,:);
% % % % %     
% % % % %     Fws = Fs(1,:)+repmat(MRF.gpe(B, qdrpGrad), size(Fs(1,:),1), 1);
% % % % %     FwI(:,:,i) = interp1(B, Fws, tempBs);
% % % % % 
% % % % % end
% % % % % 
% % % % %  figure(54)
% % % % %  for i=1:length(amp36)
% % % % %      plot(tempBs,FwI(:,:,i)); hold on
% % % % %  end
% % % % %  hold off
% % % % % 
% % % % %  % come back to this!!
% % % % % figure(155)
% % % % % striplinePlot(tempBs,(reshape(FwI,size(FwI,2), size(FwI,3))))
% % % % % % % % imagesc(amp36, Bs, FwI)
% % % % % % % % imagesc(FwI)

%% Lots of barrier heights!

% % % ladno = 3;

% rabis = [amp30 amp36(2) amp42]';
% [ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, rabis, 'iterations', 4, 'qdrpGrad', qdrpGrad);
% F = MRF.sortEnergies2(B, MRF.ladder(RFs, ladno, F)); % remove this!
% 
% Fs = zeros(size(F,1), size(F,2), length(amp36)); 

% % % for i = 1:length(amp36)
% % %     rabis = [amp30 amp36(i) amp42]';
% % %     
% % %     [ Fs, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, rabis, 'iterations', 4, 'qdrpGrad', qdrpGrad);
% % %     Fs = MRF.sortEnergies(B, MRF.ladder(RFs, ladno, Fs));
% % %     
% % %     FsTest = Fs(1,:);
% % %     
% % %     figure(100)
% % %     plot(B, i+FsTest+repmat(MRF.gpe(B, qdrpGrad), size(FsTest,1), 1),'-')
% % %     hold on
% % %     disp(i)
% % % end
% % %  hold off


 
% lad = MRF.ladder(RFs, 3, F);
% figure(1)
% plot(B, lad + repmat(MRF.gpe(B, qdrpGrad), size(lad,1), 1), '-');

% Ftest = F(1,:,1);
% 
% figure(2)
% plot(B, Ftest+repmat(MRF.gpe(B, qdrpGrad), size(Ftest,1), 1))


end
