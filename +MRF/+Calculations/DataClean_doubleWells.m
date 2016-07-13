function [ cys, navgcys, avgcys, bhs ] = DataClean_doubleWells( shots, absImgF, roi, slicex, slicey)
%DataClean Organise shell image data ready to plot
%   Creating histograms from absorption images; average repeated data and
%   sort by barrier height
%   Need to have defined absImgF


cy = sum(roi(absImgF(shots,1), slicex, slicey), 2);

% Perform for all images
cys = zeros(length(cy),size(absImgF(shots,1), 1));
for i=1:size(shots, 1)
   cys(:,i) = sum(roi(absImgF(shots,i), slicex, slicey), 2); 
end

% Sort according to barrier height.
[~,j] = sort([shots{:,1}]);
% waterfall(cys(:,j)')

% normalise total atom number in each shot.
ncys = cys ./ repmat(sum(cys, 1), size(cys, 1), 1);

% average runs with the same barrier height
exptBhs = [shots{:,1}];
bhs = unique(exptBhs);
matches = @(barrier, val) abs(barrier - val) < 1e-3; %fp precision

avgcys = zeros(size(cys, 1), length(bhs));
for i=1:length(bhs)
    barrier = bhs(i);
    t = zeros(size(cys, 1), 1);
    n = 0;
    for j=1:size(shots, 1)
        if matches(barrier, exptBhs(j))
           t = t + ncys(:,j); 
           n = n + 1;
        end
    end
    if n > 0
        t = t / n;
    end
    avgcys(:, i) = t;
end

% renormalise averages
navgcys = avgcys ./ repmat(sum(avgcys, 1), size(avgcys, 1), 1);
% waterfall(-log(navgcys'));


end

