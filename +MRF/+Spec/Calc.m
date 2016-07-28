function [ spec, debug ] = Calc( Bs, RFs, Rabi, varargin )
%CALCSPEC Calculates the spectroscopy transitions for the given MRF
%system. Outputs these transitions, plus additional diagnostic information
%that should be used to verify the calculation was correct.
% Take care: this function could give incorrect results if some of its
% assumptions falter. For example, it assumes that the minimum where the
% spec should be taken from is located within the given Zeeman splitting
% range.
% Syntax: Calc(Bs, RFs, Rabi, ... )
%  Bs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)
%
% The following parameters may also be described:
%  ladderN: length of dressed state ladder for calculating transitions.
%  iterations: number of iterations to perform meshing of eigenstates
%  qdrpGrad: quadrupole gradient to include gravity by mapping zeeman
%            splittings to height

p = inputParser;
addRequired(p,'Bs',@isnumeric);
addRequired(p,'RFs',@isnumeric);
addRequired(p,'BRFs',@isnumeric);
addParameter(p,'ladderN',5,@isnumeric);
addParameter(p,'iterations',6,@isnumeric);
addParameter(p,'qdrpGrad',0,@isnumeric);
addParameter(p,'minMethod','min',@ischar);
addParameter(p,'minStart',mean(RFs),@isnumeric);
addParameter(p,'trappedManifoldMethod','secondDeriv',@ischar);

parse(p,Bs,RFs, Rabi, varargin{:});

% Calculate the dressed energy levels and create ladder structure
[ F, B ] = MRF.MeshedQuasiEnergies(Bs, RFs, Rabi, 'iterations', p.Results.iterations, 'qdrpGrad', p.Results.qdrpGrad);
F = MRF.sortEnergies(B,MRF.ladder(RFs, 10, F));
lad = MRF.ladder(RFs, p.Results.ladderN, F);

% Select the trapped manifold. To do this, select the manifold which has a
% predominantly positive second derivative.
switch p.Results.trappedManifoldMethod
    case 'secondDeriv'
        dB = B(2:end) - B(1:end-1);
        d  = (F(:,2:end) - F(:, 1:end-1))./repmat(dB, size(F,1), 1);
        dB2 = (dB(1:end-1) + dB(2:end))/2;
        d2 = (d(:,2:end) - d(:,1:end-1))./repmat(dB2, size(d,1), 1);
        
        [~,ti] = max(mean(d2, 2));                  % trapped
        [~,uti] = min(mean(d2, 2));                 % untrapped
        ind = 1:3;
        fli = ind(~(ind == ti) & ~(ind == uti));    % flat
        
    %case 'minimaNear'
        
        % try to pick the level which has a minima nearest to our starting
        % point?
        
        
end

% Select the different manifolds. Identify the middle most trapped manifold
trapped = lad(ti:3:end, :);
flat = lad(fli:3:end, :);
antitrapped = lad(uti:3:end, :);
fav = ceil(size(trapped, 1)/2);

% Locate the potential minimum. Accuracy can be increased through iteration
% number for earlier meshing sequence.
mT = trapped(1,:);
if p.Results.qdrpGrad > 0.1
    mT = mT + MRF.gpe(B, p.Results.qdrpGrad);
end

switch p.Results.minMethod
    case 'min'
        [~,mini] = min(mT);
    case 'settle'
        mini = MRF.Spec.settleFromPoint(B, mT, p.Results.minStart);
end

spec = (flat(:,mini) - trapped(fav, mini));
spec = lad(:,mini) - trapped(fav,mini);
spec = sort(abs(spec));

debug = struct;
debug.trapped = trapped(1,:);
debug.B       = B;
debug.all = sort(abs(lad(:,mini) - trapped(fav,mini)));
debug.min = mini;

end