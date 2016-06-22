function [ lad ] = sortEnergies( B, lad, order )
%SORTENERGIES Sorts and joins the eigenstates together to remove artificial
%avoided crossings due to eg the phase folding over.
%
% To be clear:
%  Overlap: overlap of uncoupled eigen energies
%  Avoided crossing: avoided crossing between coupled eigen energies
%
% We need to be able to distinguish the two!

% F is a ladder of eigenstates of the propagator. We first identify the
% location of any overlaps, finding any points where the eigenstates come
% too closely together.

if nargin < 3
    order = 3;

thresh = 0.01;

d = (diff(lad, order, 1));
is = {};
js = {};
for j=1:size(d, 1)
    mask = abs(d(j,:)) < thresh;
    overlaps = find(mask);
    
    % overlaps may have consecutive indices if a number of points fall
    % between the threshold. This will flip-flop the eigenstates in
    % numbering and would be incorrect. We need to instead determine the
    % point at which they are nearest to each other.
    if ~isempty(overlaps)
        
        strt = [0 find(diff(overlaps) > 1)]+1;
        fnsh = [find(diff(overlaps) > 1) length(overlaps) ];
        blocks = {};
        for ind=1:length(strt)
            blocks{end+1} = overlaps(strt(ind):fnsh(ind));
        end
        
        % Within such a continous block:
        overlaps = [];
        for ind=1:length(blocks)
            [~,smallestSep] = min(abs(d(j,blocks{ind})));
            overlaps(end+1) = smallestSep + min(blocks{ind}) - 1;
        end
    end
    
    is{end+1} = overlaps;
    js{end+1} = [j j+order];
end

% We now have the locations of an intercept. We have the index of B field
% they occur immediately after and also the levels we need to swap to
% connect the eigenstates correctly.

for k=1:length(is)
    i = is{k};
    j = js{k};
    if ~isempty(i)
        for k2 = 1:length(i);
            % for each swap that needs to be made, swap them.
            temp = lad(j(1), :);
            lad(j(1), i(k2):end) = lad(j(2), i(k2):end);
            lad(j(2), i(k2):end) = temp(1,i(k2):end);
        end
    end
end

% Return only the states in the middle - ones at the very end could well be
% wrong as the states they intersect with wont have been swapped if they
% were truncated out.
middle = round(size(lad, 1)/2);
lad = lad(middle+(-1:1:1), :);

if order > 1
    lad = MRF.sortEnergies(B, lad, order-1);
end

end

