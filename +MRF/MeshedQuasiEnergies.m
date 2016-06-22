function [ eigF, Bs ] = MeshedQuasiEnergies( Bs, RFs, Rabi, varargin )
%MESHEDQUASIENERGIES Calculates the quasi energies, meshing finely over any
%points of interest (by default: stationary points)

p = inputParser;
   addRequired(p,'Bs',@isnumeric);
   addRequired(p,'RFs',@isnumeric);
   addRequired(p,'BRFs',@isnumeric);
   addParameter(p,'iterations',3,@isnumeric);
   addParameter(p,'qdrpGrad',0,@isnumeric);
   
parse(p,Bs,RFs, Rabi, varargin{:});

eigF = MRF.GetQuasiEnergies(Bs, RFs, Rabi);

iterations = p.Results.iterations;
deltaB = Bs(2)-Bs(1);

for i=1:iterations
deltaB = deltaB / 2;

cB = Util.Refine(Bs, eigF(1,:));

if p.Results.qdrpGrad > 0.1
    % if Quadrupole gradient is specified then also mesh according to gravity.
    % convert Zeeman splitting into microns, then calculate gpe.
    gp = MRF.gpe(Bs, p.Results.qdrpGrad);
    
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

