function [ minPos ] = locateWells( Bs, FworkingI, cutoff)
%locateWells Find the B-field location and corresponding energy of the trap
% minima
%   ** need to fix this so it doesn't fall apart with different numbers
%   of wells!


%% first get an estimate of the locations of the trap minima

% estimate first and second derivatives
d1 = diff(FworkingI);
d2 = diff(FworkingI,2);

% corresponding fields
Bs1 = (Bs(1:length(Bs)-1)+Bs(2:length(Bs)))/2;
Bs2 = (Bs1(1:length(Bs1)-1)+Bs1(2:length(Bs1)))/2;

% figure(3)
% % plot(Bs, 0.005*FworkingI, '.'); hold on;
% plot(Bs,zeros(length(Bs))); hold on;
% plot(Bs1, d1, '.'); hold on;
% plot(Bs2, d2, '.'); hold off

% figure(4)
% plot(Bs, FworkingI, '.');

% get rough position of zero crossings 
maxs = [];
mins = [];
stats = [];
for i = 1:length(Bs1)-1
    if d1(i)*d1(i+1)<0 && d2(i)<0
        maxs = [maxs i];
    elseif d1(i)*d1(i+1)<0 && d2(i)>0
        mins = [mins i];
    elseif d1(i)*d1(i+1)<0
        stats = [stats i];
    end
end

% maxs 
% mins
% stats
% Bs(maxs)
% Bs(mins)


% now want to search for actual minimum near these points - take 1/8
% distance between them to work with - this is quite phenomenological
% though since just trying to keep away from sharp turnover... needs to be
% more robust!

tempMin = [Bs(mins)' FworkingI(mins)'];

% % % NOTE - this definitely needs tidying!!! need to check how applicable it
% % % is... Assume for now that do only get 2 minima (since do here...)
% % if size(tempMin,1)~=2
% %     disp('Caution! Not detecting 2 minima.')
% % end
maxs;
mins;

if size(tempMin, 1) == 2
    minRange = abs(diff(mins))/cutoff % looking at divide between minima
    disp('2 minima detected')
elseif size(tempMin, 1) == 1
%     abs(maxs-mins)
    minRange = min(abs(maxs-mins));
    % choose range according to proximity to maxima - could take all the 
    % way to the maxima if wanted to map potential!
    disp('1 minimum detected')
else
    disp('No potential minimum! We have serious problems here!')
end


minRanges = zeros(round(minRange)*2+1, 2, size(tempMin,1));
minPos = zeros(size(tempMin,1), 2);

for i=1:size(tempMin,1)
    minRanges(:,:,i) = ...
        [(Bs(mins(i)-round(minRange) : mins(i)+round(minRange)))' ...
        (FworkingI(mins(i)-round(minRange) : mins(i)+round(minRange)))'];
    
    slice = minRanges(:,:,i);
    [~, index] = sort(slice(:,2));
    sorted = slice(index,:);
    minPos(i,:) = sorted(1,:);
end


end

