function [ F ] = sortEnergies2( B, lad )
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

% This second method picks three states near the centre of the ladder. It
% attempts to join the seams by finding the eigenstate that matches its
% 'best guess' for what the next value should be.

middle = round(size(lad, 1)/2);

F = lad(middle+(-1:1:1), :);
dB = diff(B);
dE = diff(lad, 1, 2);

bests = [];
guesses = [];
% iterate over three states near the centre of the ladder
for i=(1:3)
    mi = middle+i-2;
    guesses(end+1) = F(i,1);
    % iterate by running over the eigenstate values
    for j=3:size(lad, 2)-1
       
        % guess what next value should be based on previous
        guess = lad(mi, j-1) + dE(mi,j-2)./dB(j-2).*dB(j-1);
        
        % find the eigenstate that best matches the guess
        [~,mi] = min(abs(guess - lad(:, j)));
        
        % write this new value into our eigenstate
        F(i, j) = lad(mi, j);
        bests(end+1) = mi;
        guesses(end+1) = guess;
    end
    
    F(i, end) = lad(mi, end);
    guesses(end+1) = guess;
end

end
