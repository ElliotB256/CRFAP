function [ E, Ev, barestates ] = dressedStateBasis( Bs, RF, N, g )
%DRESSEDSTATEBASIS Finds the dressed state eigenvalues and vectors

ip = inputParser;

addRequired(ip, 'Bs'); % Zeeman splitting in MHz
addRequired(ip, 'RF'); % dressing RF. Can be multi component.
addRequired(ip, 'N');  % number of photon modes. Can be multi component.
addRequired(ip, 'g');  % interaction strength. Can be multi component.

parse(ip, Bs, RF, N, g);

H_a = diag( Bs * [ -1 0 1 ]' );
H_RF1 = diag( RF(1) * (1:N(1)) );

if length(RF) > 1
    H_RF2 = diag( RF(2) * (1:N(2)));
else
    H_RF2 = 1;
end

% identity matrices:
I_a      = eye(3);
I_RF1    = eye(N(1));

if length(RF) > 1
    I_RF2    = eye(N(2));
else
    I_RF2 = 1;
end

% define basis: occupation numbers for the bare states are easy to use.
if length(N) > 1
    barestates = labelStates( [ 3 N(1) N(2) ] );
else
    barestates = labelStates( [ 3 N(1) 1 ] );
end

merge = @(a,b,c) kron(I_RF2, kron(I_RF1, a)) + kron(I_RF2, kron(b, I_a)) + kron(c, kron(I_RF1, I_a));

% Calculate the Hamiltonian of the atomic energy plus the field energies.
H0 = merge(H_a, H_RF1, H_RF2);

H_int = zeros(size(barestates, 1));
% Construct interaction matrix in a more general way using labelled states
for i=1:size(H_int,1)
    for j=1:size(H_int,2)
        
        a = barestates(i, :);
        b = barestates(j, :);
        
        % both are raised with no effect on probe N:
        if a(1) == b(1) + 1 && a(2) == b(2) + 1 && a(3) == b(3)
            H_int(i,j) = 1;
        end
        
        % both are lowered with no effect on probe N:
        if a(1) == b(1) - 1 && a(2) == b(2) - 1 && a(3) == b(3)
            H_int(i,j) = 1;
        end
    end
end
H_int = H_int * g(1);

% We now can formulate the complete Hamiltonian
H = H0 + H_int;
%imagesc(H);

% Calculate the eigenstates of this Hamiltonian
[E,Ev] = eig(H);

% Define transformation from 'bare' basis to dressed state basis
% toDressedBasis = @(U) E \ U * E;

end

