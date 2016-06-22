function [ spec, debug ] = Calc( Bs, RFs, Rabi, varargin )
%CALCSPEC Calculates the spectroscopy transitions for the given MRF
%system. Outputs these transitions, plus additional diagnostic information
%that should be used to verify the calculation was correct. 
% Take care: this function could give incorrect results if some of its
% assumptions falter. For example, it assumes that the minimum where the
% spec should be taken from is located within the given Zeeman splitting
% range.

p = inputParser;
   addRequired(p,'Bs',@isnumeric);
   addRequired(p,'RFs',@isnumeric);
   addRequired(p,'BRFs',@isnumeric);
   addParameter(p,'ladderN',5,@isnumeric);
   addParameter(p,'iterations',6,@isnumeric);
   
parse(p,Bs,RFs, Rabi, varargin{:});

% Calculate the dressed energy levels and create ladder structure
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', p.Results.iterations);
lad = MRF.ladder(RFs, p.Results.ladderN, F);

% Select the trapped manifold. To do this, select the manifold which has a
% predominantly positive second derivative.

d2 = diff(diff([F F(:,end)]'))'./repmat(diff(B).^2, size(F,1), 1);
[~,ti] = max(mean(d2, 2));                  % trapped
[~,uti] = min(mean(d2, 2));                 % untrapped
ind = 1:3;
fli = ind(~(ind == ti) & ~(ind == uti));    % flat

% Select the different manifolds. Identify the middle most trapped manifold
trapped = lad(ti:3:end, :);
flat = lad(fli:3:end, :);
antitrapped = lad(uti:3:end, :);
fav = ceil(size(trapped, 1)/2);

% Plot the energy level structure. Identify the trapped manifold; important
% to check that this is the correct one!
%figure;
%plot(B, trapped, 'r'); hold on; plot(B, flat, 'Color', [0.5 0.5 0.5]); plot(B, antitrapped, 'Color', [0 0 0]); plot(B, trapped(fav, :), 'r', 'LineWidth', 2); hold off
%plot(B, F(ti,:), 'r'); hold on; plot(B, F(fli,:), 'Color', [0.5 0.5 0.5]); plot(B, F(uti,:), 'Color', [0 0 0]); hold off;

% Locate the potential minimum. Accuracy can be increased through iteration
% number for earlier meshing sequence.
[~,mini] = min(trapped(1,:));

spec = (flat(:,mini) - trapped(fav, mini));
%spec = lad(:,mini) - trapped(fav,mini);
spec = sort(abs(spec));

debug = struct;
debug.trapped = trapped(1,:);
debug.B       = B;
debug.all = sort(abs(lad(:,mini) - trapped(fav,mini)));

end