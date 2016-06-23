function [ F ] = sortEnergies( B, lad )
%SORTENERGIES Sorts and joins the eigenstates together to remove artificial
%avoided crossings due to eg the phase folding over or overlap of states.
% This method uses a trajectory based method, where we travel along each
% state and find the eigenvalues that best match a guess for the next
% value based on previous curvature of the current state.
%
% Syntax: sortEnergies(B, lad)
%  B: energy splitting of the undressed Zeeman states in MHz.
%  lad: ladder of dressed states. This ladder should be great enough that
%       truncation does not give erroneous results (try increasing size if 
%       results look incorrect!).



% To be clear:
%  Overlap: overlap of uncoupled eigen energies
%  Avoided crossing: avoided crossing between coupled eigen energies
% 
% F is a ladder of eigenstates of the propagator. We first identify the
% location of any overlaps, finding any points where the eigenstates come
% too closely together.

% This second method picks three states near the centre of the ladder. It
% attempts to join the seams by finding the eigenstate that matches its
% 'best guess' for what the next value should be.

middle = round(size(lad, 1)/2);

F = lad(middle+(-1:1:1), :);
%dE = lad(:,2:end) - lad(:, 1:end-1);
%dB = B(2:end) - B(1:end-1);

guesses = F;
% store of best matched indices
best = repmat((1:3)', 1, size(F, 2));

% iterate over three states near the centre of the ladder
for i=(1:3)
    mi = middle+i-2;
    guesses(i, 1:2) = F(i,1:2);
    best(i, 1:2) = mi;
    F(i, 1:2) = lad(mi, 1:2);

    % iterate by running over the eigenstate values
    for j=2:size(lad, 2)-2
       
        % Need to use present values for trajectory: unfortunately can't
        % precompute diffs easily.
        
        % guess what next value should be based on previous.        
        pi = best(i,j-1);
        guess = lad(mi, j) + (B(j+1) - B(j)) ./ (B(j) - B(j-1)) .* (lad(mi, j) - lad(pi,j-1));
        
        % find the eigenstate that best matches the guess
        [~,mi] = min(abs(guess - lad(:, j+1)));
        
        % write this new value into our eigenstate
        F(i, j+1) = lad(mi, j+1);
        best(i, j+1) = mi;
        guesses(i,j+1) = guess;
    end

    F(i, end) = lad(mi, end);
    best(i, end) = mi;
end

end
