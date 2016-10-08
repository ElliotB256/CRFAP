function [ F ] = sortEnergies( Zs, lad, varargin )
%SORTENERGIES Sorts and joins the eigenstates together to remove artificial
%avoided crossings due to eg the phase folding over or overlap of states.
% This method uses a trajectory based method, where we travel along each
% state and find the eigenvalues that best match a guess for the next
% value based on previous curvature of the current state.
%
% Syntax: sortEnergies(B, lad)
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  lad: ladder of dressed states. This ladder should be great enough that
%       truncation does not give erroneous results (try increasing size if 
%       results look incorrect!).

p = inputParser;
addRequired(p,'Zs',@isnumeric);
addRequired(p,'lad',@isnumeric);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));
parse(p, Zs, lad, varargin{:});


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

switch p.Results.F
    case 1
        spaceSize = 3;
        md = -1:1:1;
    case 2
        spaceSize = 5;
        md = -2:1:2;
end

F = lad(middle+md, :);
guesses = F;
% store of best matched indices
best = repmat((1:spaceSize)', 1, size(F, 2));

% iterate over three states near the centre of the ladder
for i=(1:spaceSize)
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
        guess = lad(mi, j) + (Zs(j+1) - Zs(j)) ./ (Zs(j) - Zs(j-1)) .* (lad(mi, j) - lad(pi,j-1));
        
        % guess also what DERIVATIVE should be next based on the previous
        % values
        
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
