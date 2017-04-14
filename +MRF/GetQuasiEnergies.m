function [ eigF2, eigV ] = GetQuasiEnergies( Zs, RF, gFuBB, varargin )
%GETQUASIENERGIES Calculates quasi-energies for the given multi-RF
%field over the specified range of zeeman energy splittings.
% Syntax: GetQuasiEnergies( Zs, RF, Rabi)
%  Zs: energy splitting of the undressed Zeeman states in MHz.
%  RF: vector of dressing RFs (MHz)
%  gFuBB: vector of dressing RF amplitudes (MHz)

p = inputParser;
addRequired(p,'Zs',@isnumeric);
addRequired(p,'RF',@isnumeric);
addRequired(p,'gFuBB',@isnumeric);
addParameter(p,'phase',0,@isnumeric);
addParameter(p,'F', 1, @(x) any(ismember(x,[1 2])));
addParameter(p, 'theta', 0, @(x) all(x >= 0 & x <= pi) && size(x, 2) == 1);
addParameter(p, 'polarisation', 'circ', @(x) ismember(x, {'circ', 'lin'}));
addParameter(p, 'parallel', 1); % Try to evaluate in parallel
parse(p,Zs, RF, gFuBB, varargin{:});

switch p.Results.F
    case 1
        hs = 3;
    case 2
        hs = 5;
end

% Only calculate eigen values if we request them.
calcEigV = nargout > 1;

eigF2 = zeros(hs, length(Zs));
periodicity = 2*pi/MRF.GetFundamental(RF);

if calcEigV
    eigV = zeros(hs, hs, length(Zs));
end

% Try to perform the calculations in parallel where possible.
parfailed = 0;
if p.Results.parallel
    try
        parfor i=1:length(Zs)
            B = Zs(i);
            if ~calcEigV
                eigF2(:,i) = compute(B, RF, gFuBB, p, periodicity);
            else
                [ eigF2(:,i), eigV(:,:,i) ] = compute(B, RF, gFuBB, p, periodicity);
            end
        end
    catch e
        
        if strcmp(e.identifier, 'parallel:lang:pool:UnexpectedParforFailure')
            fprintf('Unable to checkout parallel toolbox license. Falling back to non-parallel computation.\n')
            parfailed = 1;
        end
        
    end
end

% If parallel execution is not possible, fall back to standard.
if ~p.Results.parallel || parfailed
    for i=1:length(Zs)
        B = Zs(i);
        if ~calcEigV
            eigF2(:,i) = compute(B, RF, gFuBB, p, periodicity);
        else
            [ eigF2(:,i), eigV(:,:,i) ] = compute(B, RF, gFuBB, p, periodicity);
        end
    end
end

end