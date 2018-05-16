function Sample(instance)
%SAMPLE Calculate eigenenergies at each B.

lambda = instance.InitialLambda;
lambdaTol = min(diff(lambda));

[B,theta,gamma] = instance.lambdaToBTG(lambda);

if (instance.Verbose)
    fprintf('Calculating initial energies at %d field points.\n', length(B))
end
% Calculate first set of eigenvalues at initial points
[eigE,eigV] = instance.APCalculator.GetDressedEnergies(B, theta, gamma);

for i=1:instance.MeshIterations
    
    cl = Util.Refine(lambda, eigE(1,:));
    
    if instance.QuadGrad > 0.0001 && instance.Theta ~= pi/2
        % if Quadrupole gradient is specified then also mesh according to gravity.
        % convert Zeeman splitting into microns, then calculate gpe.
        
        % Get z of each point and apply GPE
        coords = instance.TransformThetaGamma2XYZ(instance.lambdaToBTG(lambda), instance.Theta, instance.Gamma);
        gp = Util.gpe(coords.z, instance.APCalculator.Atom.Mass);
        % Modify eigen energies to account for energy shift
        modEig = eigE + repmat(gp, instance.APCalculator.HilbertSpaceSize, 1);
        
        % Collect unique values suggested from eigenvalues
        gcl = [];
        for j=1:size(eigE, 1)
            gcl = [gcl Util.Refine(lambda, modEig(j,:))];
        end
        gcl = uniquetol(gcl, lambdaTol/20);
        
        % Remove duplicate elements
        cl = Util.UniquePick( lambda, [ cl gcl ]);
    end
    
    % Calculate energies of these new points and add to the list.
    lambda2 = cl(:)';
    if (instance.Verbose)
        fprintf('Iteration %d of %d: Calculating meshed energies at %d field points.\n', i, instance.MeshIterations, length(lambda2))
    end
    
    [B,theta,gamma] = instance.lambdaToBTG(lambda2);
    [eigE2,eigV2] = instance.APCalculator.GetDressedEnergies(B, theta, gamma);
    eigE = [eigE eigE2];
    eigV = cat(3, eigV, eigV2);
    lambda = [lambda lambda2];
    
    % sort results by B
    [~,j] = sort(lambda);
    lambda = lambda(j);
    eigE = eigE(:,j);
    eigV = eigV(:,:,j);
    
end

instance.UnsortedEigenenergies = eigE;

% Sort eigenvalues using similarity between eigenvectors.
if (instance.Sort)
    
    if (instance.Verbose)
        fprintf('Sorting eigenvectors at %d field points.\n', length(B))
    end
    
    eigE2 = zeros(size(eigE));
    eigE2(:,1) = eigE(:,1);
    eigV2 = zeros(size(eigV));
    eigV2(:,:,1) = eigV(:,:,1);
    currIndices = 1:size(eigE);
    
    % Sweep across the spatial coordinate. At each point, associate each
    % eigenvector from the previous coordinate with the eigenvector that
    % has the most overlap.
    %
    % Note that we also have to look at the negative of the eigenstates, as
    % the sign can change and it is still really the same eigenstate.
    %
    % Output of sorting will be stored in eigE2.
    
    for i=2:length(lambda)
        for j=1:size(eigV, 1)
            % determine the eigenstate which is most similar to the previous.
            
            prevState = squeeze(eigV(:,currIndices(j),i-1));
            candidates = eigV(:,:,i);
            difference = @(a,b) sum(conj(a).*b, 1);
            scores = [ difference(repmat(prevState, 1, size(eigV,1)), candidates) difference(repmat(-prevState, 1, size(eigV,1)), candidates)];
            %             scores = difference(repmat(prevState, 1, size(eigV,1)), candidates);
            [~,next] = max((scores));
            if next > size(eigV,1)
                next = next-size(eigV,1);
            end
            currIndices(j) = next;
            
            eigV2(:,j,i) = eigV(:,next,i);
            
            % Having determined the eigenstate that is most similar, find the
            % eigenenergy of that state +- periodicity that is nearest to
            % previous energy.
            prevE = eigE2(j, i-1);
            nextE = eigE(next, i);
            
            fundamental = Floquet.GetFundamental(instance.APCalculator.RF);
            candidates = nextE + (-100:100)*fundamental;
            [~,k] = min(abs(candidates - prevE));
            nextE = candidates(k);
            
            eigE2(j, i) = nextE;
        end
    end
    
    % Use sorted values
    eigE = eigE2;
    eigV = eigV2;
    
end

%Store results in sampler object.
instance.mLambda = lambda;
instance.mE = eigE;
instance.mEV = eigV;
instance.Dirty = 0;

end