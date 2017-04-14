%% mF Decomposition
% Calculate the mF composition of eigenvectors of the propagator T, the
% application of which evolves the wavefunction for one period of time. As
% this operator T commutes with the Hamiltonian, these eigenvectors are
% also the eigenvectors of the dressed Hamiltonian.
RF = [ 3 3.6 4.2 ]';
Rabi = [ 0.1 0.1 0.1 ]';
qdrpGrad = 100;

Bs=linspace(2.5, 4.7, 50);
[ eigF, eigV ] = MRF.GetQuasiEnergies(Bs, RF, Rabi, 'parallel', 0);
%%
% Let's try and sort the eigenstates and stitch together the linked ones.
% Rather than using the usual method of adiabatically joining the states
% together, let's take an approach based on the eigenvectors. We choose the
% next value of a given eigenstate by choosing the eigenvector that is most
% similar to the previous value of it's eigenvector.
difference = @(a,b) sum((a-b).*conj(a-b), 1);

eigF2 = zeros(size(eigF));
eigF2(:,1) = eigF(:,1);
eigV2 = zeros(size(eigV));
eigV2(:,:,1) = eigV(:,:,1);
currIndices = [ 1 2 3 ];

% i indexes the horizontal position we are in.
for i=2:length(Bs)
    % j indexes the eigenstates.
    for j=1:3
        % determine the eigenstate which is most similar to the previous.
        % (we also have to look at the negative of the eigenstates, as the
        % sign can change and it is still really the same eigenstate)
        prevState = squeeze(eigV(:,currIndices(j),i-1));
        candidates = eigV(:,:,i);
        scores = [ difference(repmat(prevState, 1, 3), candidates) difference(repmat(-prevState, 1, 3), candidates)];
        [~,next] = min(scores);
        if next > 3
            next = next-3;
        end
        lastIndex = currIndices(j);
        currIndices(j) = next;
        
        eigV2(:,j,i) = eigV(:,next,i);
        
        % Having determined the eigenstate that is most similar, find the
        % eigenenergy of that state +- periodicity that is nearest to
        % previous energy.
        prevE = eigF2(j, i-1);
        nextE = eigF(next, i);
        
        periodicity = 0.6;
        candidates = nextE + (-100:100)*periodicity;
        [~,k] = min(abs(candidates - prevE));
        nextE = candidates(k);
        
        eigF2(j, i) = nextE;
    end
end


% Plot the colored line in Matlab. We have to use a patch object as a
% straight plot won't support per-vertex color data.
for i=1:3
    x = [ Bs NaN ]; %NaN are used to prevent the patch looping back
    y = [eigF2(i,:) NaN ];
    c = squeeze(abs(eigV2(:,i,:)))'; c = [ c; NaN NaN NaN ];
    p = patch(x,y,x, 'CDataMapping', 'direct', 'FaceColor', 'interp', 'EdgeColor', 'interp', 'FaceAlpha', 0);
    p.FaceVertexCData = c; p.LineWidth = 2; hold on;
end
hold off;

xlabel('Zeeman splitting (MHz)', 'FontSize', 14);
ylabel('V/h (MHz)', 'FontSize', 14);

set(gcf, 'Color', 'w');
box on;