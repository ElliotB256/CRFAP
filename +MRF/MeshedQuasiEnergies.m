function [ eigF, Bs ] = MeshedQuasiEnergies( Bs, RFs, Rabi, varargin )
%MESHEDQUASIENERGIES Calculates quasi-energies for the given multi-RF
%field over the specified range of zeeman energy splittings. The values Bs
%specify the points where the potential is evaluated initially. The result
%is improved over a given number of iterations by meshing over interesting
%features.
% Syntax: GetQuasiEnergies( Bs, RFs, Rabi, ...)
%  Bs: energy splitting of the undressed Zeeman states in MHz.
%  RFs: vector of dressing RFs (MHz)
%  Rabi: vector of dressing RF Rabi frequencies (MHz)
% 
% The following parameters may also be described:
%  iterations: number of times to more finely mesh the results.
%  qdrpGrad: incorporates gravity using the specified qdrpGrad to convert
%            frequency to space

p = inputParser;
   addRequired(p,'Bs',@isnumeric);
   addRequired(p,'RFs',@isnumeric);
   addRequired(p,'BRFs',@isnumeric);
   addParameter(p,'iterations',3,@isnumeric);
   addParameter(p,'qdrpGrad',0,@isnumeric);
   addParameter(p,'gF',Constants.gF,@isnumeric);
   addParameter(p,'mass',87,@isnumeric);
   
parse(p,Bs,RFs, Rabi, varargin{:});

eigF = MRF.GetQuasiEnergies(Bs, RFs, Rabi);

iterations = p.Results.iterations;
deltaB = Bs(2)-Bs(1);

for i=1:iterations
deltaB = deltaB / 2;

cB = Util.Refine(Bs, eigF(1,:));

if p.Results.qdrpGrad > 0.0001
    % if Quadrupole gradient is specified then also mesh according to gravity.
    % convert Zeeman splitting into microns, then calculate gpe.
    gp = MRF.gpe(Bs, p.Results.qdrpGrad, 'gF', p.Results.gF, 'mass', p.Results.mass);
    
    % Modify eigen energies to account for energy shift
    modEig = eigF + repmat(gp, 3, 1);
    
    % Now collect unique values suggested from the first and third
    % eigenvalues
    gcB = Util.Refine(Bs, modEig(1,:));
    gcB2 = Util.Refine(Bs, modEig(3,:));
    
    % Remove duplicate elements
    cB = Util.UniquePick( Bs, [cB gcB gcB2 ]);
end

% Calculate energies of these new points and add to the list.
Bs2 = cB(:)';
eigF2 = MRF.GetQuasiEnergies(Bs2, RFs, Rabi);
eigF = [eigF eigF2];
Bs = [Bs Bs2];

% sort results by B
[~,j] = sort(Bs);
Bs = Bs(j);
eigF = eigF(:,j);

end
end

