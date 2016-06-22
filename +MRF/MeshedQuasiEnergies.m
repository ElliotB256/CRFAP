function [ eigF, Bs ] = MeshedQuasiEnergies( Bs, RFs, Rabi, varargin )
%MESHEDQUASIENERGIES Calculates the quasi energies, meshing finely over any
%points of interest (by default: stationary points)

p = inputParser;
   addRequired(p,'Bs',@isnumeric);
   addRequired(p,'RFs',@isnumeric);
   addRequired(p,'BRFs',@isnumeric);
   addParameter(p,'iterations',3,@isnumeric);
   
parse(p,Bs,RFs, Rabi, varargin{:});

eigF = MRF.GetQuasiEnergies(Bs, RFs, Rabi);

iterations = p.Results.iterations;
deltaB = Bs(2)-Bs(1);

for i=1:iterations
deltaB = deltaB / 2;
% locate interesting points to mesh around. Here we find stationary points.
%d = diff(eigF, 1, 2)';
%d2 = diff(sign(d), 1, 1);
%j = find(abs(d2(:,1)) > 0.5);
%j2 = min(j+1, length(Bs));



% mesh around the stationary points
% mesh = [-deltaB deltaB];
%cB = unique([Bs(j) Bs(j2)]);
% cB = unique(mean([Bs(j);Bs(j2)], 1));
%Bs2 = repmat(cB, size(mesh, 2), 1) + repmat(mesh', 1, size(cB, 2));
%Bs2 = Bs2(:)';
% Bs2 = cB(:)';

cB = Util.Refine(Bs, eigF(1,:));

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

