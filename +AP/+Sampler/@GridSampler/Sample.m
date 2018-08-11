function Sample(sampler)
%SAMPLE Calculates dressed eigenenergies at each point in the sampler.

x = sampler.X(:);
y = sampler.Y(:);
z = sampler.Z(:);

sampler.B = (sampler.QuadGrad * 1e-4) .* (x.^2 + y.^2 + 4 * z.^2).^0.5;
% First, calculate distribution of theta and gamma used for stated x,y,z
% points
sampler.Gamma = acos(x ./ (x.^2 + y.^2).^0.5);
sampler.Theta = acos(((x.^2 + y.^2) ./ (x.^2 + y.^2 + 4 * z.^2)).^0.5);

if sampler.Verbose
    fprintf('Calculating %d eigenstates...\n', length(sampler.B))
end

[eigE,eigV] = sampler.APCalculator.GetDressedEnergies(sampler.B, sampler.Theta, sampler.Gamma);


% Meshing procedure -
%  1. Create delaunay triangulation of points.
%  2. Start at point (0,0), sort neighbouring eigenstates according to
%  those most like that at current, iterate through tree.
if sampler.Verbose
    fprintf('Meshing %d eigenstates...\n', length(sampler.B))
end

d = delaunayTriangulation([x,y,z]);

% Enumerate until all points are assigned states
states = NaN(size(eigE,1),size(d.Points,1));
sortedEigE = eigE;
sortedEigV = eigV;

% start from first element in list.
queue = 1; states(1:size(states, 1), 1) = 1:size(states,1);
while any(any(isnan(states)))
    i = queue(1); queue(1) = [];
    
    % Get all points connected to i
    mask = any(d.ConnectivityList == i, 2);
    points = d.ConnectivityList(mask,:);
    points = unique(points(:));
    
    % If points have not had states assigned, assign based on values for i
    for j=points'
        if ~isnan(states(1,j))
            continue;
        end
        
        % point 1189
        
        % determine candidate with best overlap by maximising score
        a = sortedEigV(:,:,i);
        b = eigV(:,:,j);
        
        candidates = perms(1:size(states,1))';
        scores = zeros(length(candidates), 1);
        for k=1:length(scores)
            a0 = a(:,states(:,i));
            b0 = b(:,candidates(:,k));
            scores(k) = sum(abs(diag(conj(a0')*b0)));
        end
        [~,k] = max(scores);
        states(:,j) = candidates(:,k);
        sortedEigV(:,:,j) = b(:,states(:,j));
        
        % Assign eigenenergy that is closest to original point, within periodicity of system.
        period = Floquet.GetFundamental(sampler.APCalculator.RF);
        ladder = (-100:100)*period;
        for m=1:size(sortedEigE,1)
            difference = sortedEigE(m, i) - eigE(states(m,j), j);
            
            [~,n] = min(abs(difference + ladder));
            sortedEigE(m,j) = eigE(states(m,j),j) - ladder(n);
        end
        
        % assign the neighbour to queue to check for more connections
        queue(end+1) = j;
    end
    
end

sampler.Eigenenergies = sortedEigE;
sampler.Eigenvectors = sortedEigV;

gs = sampler
mask = abs(gs.Y) < 1e-4;
E = gs.Eigenenergies(1,:); %PotentialEnergies(1,:);
tris = delaunay(gs.X(mask), gs.Z(mask));
trisurf(tris, gs.X(mask), gs.Y(mask), gs.Z(mask), E(mask), 'EdgeColor', 'w', 'LineStyle', '-', 'LineWidth', 1); shading interp; %axis equal



if sampler.Verbose
    fprintf('Complete.\n')
end

end